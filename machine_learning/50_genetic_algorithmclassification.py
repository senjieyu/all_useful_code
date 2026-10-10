# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 09:47:41 2025

@author: senji
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, LeaveOneOut, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from deap import base, tools, creator, algorithms
import random
import os

# 创建保存结果的目录
results_dir = "results"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# 读取数据
dataset = pd.read_excel("D:/公司/卢天-纪晓波25-3-17/纪晓波20250216/分类问题/特征筛选mrmr_selected.xlsx")

# 检查目标变量的分布
unique, counts = np.unique(dataset['Class-3'], return_counts=True)
print("Class distribution:", dict(zip(unique, counts)))

# 分离特征和目标变量
X = dataset.drop(columns=['Class-3'])  # 特征
Y = dataset['Class-3']  # 目标变量

# 定义 validate 函数
def validate(X, Y, Xtest=None, Ytest=None, task="train", cv=5, model=None, **model_p):
    predictions = []
    observations = []
    if task == "train":
        model_ = model(**model_p).fit(X, Y)
        predictions = model_.predict(X)
        observations = Y.copy()
    elif task == "test":
        model_ = model(**model_p).fit(X, Y)
        predictions = model_.predict(Xtest)
        observations = Ytest.copy()
    elif task == "loo":
        loo = LeaveOneOut()
        for train_index, test_index in loo.split(X):
            xtrain, ytrain = X[train_index], Y[train_index]
            xtest, ytest = X[test_index], Y[test_index]
            model_ = model(**model_p).fit(xtrain, ytrain)
            predictions += model_.predict(xtest).tolist()
            observations += ytest.tolist()
    elif task == "cv":
        cv = KFold(cv)
        for train_index, test_index in cv.split(X):
            xtrain, ytrain = X[train_index], Y[train_index]
            xtest, ytest = X[test_index], Y[test_index]
            model_ = model(**model_p).fit(xtrain, ytrain)
            predictions += model_.predict(xtest).tolist()
            observations += ytest.tolist()
    accuracy = accuracy_score(observations, predictions)
    misclassification_rate = 1 - accuracy
    return misclassification_rate, accuracy

# 定义 fitness 函数
def fitness(individual, X, Y):
    if sum(individual) <= 2:
        return 9999,
    else:
        X_ = X[:, individual]
        return validate(X=X_, Y=Y, task="train", model=model,**model_p)[0],

# 创建遗传算法工具箱
creator.create("FitnessMax", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)
toolbox = base.Toolbox()
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=10)
toolbox.register("gen_idx", random.randint, 0, X.shape[1] - 1)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.gen_idx, n=5)  # 选择 5 个特征
toolbox.register("population", tools.initRepeat, list, toolbox.individual, n=35)

# 循环 50 次
for i in range(50):
    print(f"Running iteration {i + 1}...")

    # 使用分层抽样划分数据集
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=i, stratify=Y)

    # 标准化特征
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # 设置模型和参数
    model = SVC
    model_p = {}

    # 注册 evaluate 函数
    toolbox.register("evaluate", fitness, X=X_train, Y=Y_train)

    # 运行遗传算法
    population = toolbox.population()
    halloffame = tools.HallOfFame(1)
    pop, logbook = algorithms.eaSimple(population, toolbox, cxpb=0.5, mutpb=0.7, ngen=20, halloffame=halloffame, verbose=False)

    # 获取最佳特征子集
    fmask = halloffame[0]
    filtered_X_train = X_train[:, fmask]
    filtered_X_test = X_test[:, fmask]

    # 评估模型
    train_result = validate(X=filtered_X_train, Y=Y_train, task="train", model=model,**model_p)
    test_result = validate(X=filtered_X_train, Y=Y_train, task="test", model=model,**model_p, Xtest=filtered_X_test, Ytest=Y_test)

    # 保存结果到 CSV 文件
    result_df = pd.DataFrame({
        "Iteration": [i + 1],
        "Best_Feature_Mask": [fmask],
        "Train_Misclassification_Rate": [train_result[0]],
        "Train_Accuracy": [train_result[1]],
        "Test_Misclassification_Rate": [test_result[0]],
        "Test_Accuracy": [test_result[1]]
    })

    result_file = os.path.join(results_dir, f"result_{i + 1}.csv")
    result_df.to_csv(result_file, index=False)

    print(f"Results saved to {result_file}")

print("All iterations completed!")