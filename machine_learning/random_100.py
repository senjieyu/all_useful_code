# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 17:32:36 2024

@author: senji
"""

import pandas as pd

from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
import numpy as np



# 读取训练集和测试集
df_train = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林训练集.csv')
df_test = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林测试集.csv')

# 合并两个数据集
df_combined = pd.concat([df_train, df_test], ignore_index=True)

# 提取特征列和目标值列
X = df_combined.iloc[:, 1:]  # 假设特征列从第二列开始
y = df_combined.iloc[:, 0]   # 假设目标值列为第一列


from sklearn.metrics import accuracy_score

# 初始化随机森林分类器
# rf_classifier = RandomForestClassifier(max_depth= 5, min_samples_leaf= 1, min_samples_split=3, n_estimators= 100)
rf_classifier = RandomForestClassifier(max_depth= None, min_samples_split= 4, n_estimators= 200, random_state= 325)
# rf_classifier = RandomForestClassifier(random_state=42)
# 定义划分比例
test_size = 0.2  # 测试集比例
n_splits = 100   # 重复次数

# 初始化训练集和测试集的准确率列表
train_accuracy_scores = []
test_accuracy_scores = []

# 执行100次随机划分和模型训练
for _ in range(n_splits):
    # 随机划分训练集和测试集
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size)
    
    # 训练模型
    rf_classifier.fit(X_train, y_train)
    
    # 计算训练集和测试集的准确率
    train_accuracy = accuracy_score(y_train, rf_classifier.predict(X_train))
    test_accuracy = accuracy_score(y_test, rf_classifier.predict(X_test))
    
    # 将准确率添加到列表中
    train_accuracy_scores.append(train_accuracy)
    test_accuracy_scores.append(test_accuracy)

# 计算训练集和测试集准确率的平均值
average_train_accuracy = np.mean(train_accuracy_scores)
average_test_accuracy = np.mean(test_accuracy_scores)

print("Average Train Accuracy across 100 splits:", average_train_accuracy)
print("Average Test Accuracy across 100 splits:", average_test_accuracy)



from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

# 初始化训练集和测试集的准确率、精确度、召回率和F1分数列表
train_accuracy_scores = []
test_accuracy_scores = []
precision_scores = []
recall_scores = []
f1_scores = []

# 执行100次随机划分和模型训练
for _ in range(n_splits):
    # 随机划分训练集和测试集
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size)
    
    # 训练模型
    rf_classifier.fit(X_train, y_train)
    
    # 预测
    y_pred_train = rf_classifier.predict(X_train)
    y_pred_test = rf_classifier.predict(X_test)
    
    # 计算训练集和测试集的准确率
    train_accuracy = accuracy_score(y_train, y_pred_train)
    test_accuracy = accuracy_score(y_test, y_pred_test)
    
    # 计算精确度
    precision = precision_score(y_test, y_pred_test, average='weighted')
    
    # 计算召回率
    recall = recall_score(y_test, y_pred_test, average='weighted')
    
    # 计算F1分数
    f1 = f1_score(y_test, y_pred_test, average='weighted')
    
    # 将评估指标添加到列表中
    train_accuracy_scores.append(train_accuracy)
    test_accuracy_scores.append(test_accuracy)
    precision_scores.append(precision)
    recall_scores.append(recall)
    f1_scores.append(f1)

# 计算训练集和测试集准确率、精确度、召回率和F1分数的平均值
average_train_accuracy = np.mean(train_accuracy_scores)
average_test_accuracy = np.mean(test_accuracy_scores)
average_precision = np.mean(precision_scores)
average_recall = np.mean(recall_scores)
average_f1 = np.mean(f1_scores)

print("Average Train Accuracy across 100 splits:", average_train_accuracy)
print("Average Test Accuracy across 100 splits:", average_test_accuracy)
print("Average Precision across 100 splits:", average_precision)
print("Average Recall across 100 splits:", average_recall)
print("Average F1 Score across 100 splits:", average_f1)
