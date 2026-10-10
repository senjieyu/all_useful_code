# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 14:56:20 2024

@author: ma'wei'bin
"""
from fml import read_data
from fml.sampling import random_split
from sklearn.preprocessing import StandardScaler
import random
from sklearn.metrics import r2_score, mean_squared_error, accuracy_score
from sklearn.model_selection import LeaveOneOut, KFold
from sklearn.svm import SVR
from xgboost import XGBRegressor as xgb
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor as RF
from sklearn.ensemble import GradientBoostingRegressor as gbr
from sklearn.tree import DecisionTreeRegressor as DTR
from sklearn.linear_model import Ridge

for i in range(50):
    dataset = read_data("preprocessed_data.xlsx", df=False)
    dataset.X = StandardScaler().fit_transform(dataset.X)
    
    random_state = random.randint(0, 999999)
    train, test = random_split(dataset, percent=0.2, random_state=random_state)
    
    X_train = train.X.copy()
    X_test = test.X.copy()
    Y_train = train.Y
    Y_test = test.Y
    fnames = train.Xnames.copy()
    
    # 遗传算法参数
    ngen = 100 # 遗传代数
    cxpb =  0.8 # 交叉
    mutpb = 0.05 # 感染
    population = 50 # 种群大小
    ga_task = "cv" # 搜索时验证方法，loo、cv、train
    ga_cv = 5 # 如果ga_task为cv，那此处表示cv5
    output_feature_number = random.randint(5, 16)
    model = gbr
    model_p = dict(
        # n_estimators=100,
        # learning_rate=0.1,
        # max_depth=3,
        #alpha=0.1
        )
    
    from deap import base, tools, creator, algorithms
    
    def validate(X, Y, Xtest=None, Ytest=None, task="train", cv=5, problem="reg", model=None, **model_p):
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
        if problem == "reg":
            r2 = r2_score(observations, predictions)
            rmse = mean_squared_error(observations, predictions) ** 0.5
            return rmse, r2
        elif problem == "cls":
            return accuracy_score(observations, predictions)
    
    def fitness(individual, X, Y):
        if sum(individual) <= 2:
            return 9999,
        else:
            X_ = X[:, individual]
            return validate(X=X_, Y=Y, task=ga_task, cv=ga_cv, model=model, **model_p)[0],
    
    def gen(max_features, min_features=0):
        return random.randint(min_features, max_features-1)
    
    creator.create("FitnessMax", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)
    toolbox = base.Toolbox()
    toolbox.register("evaluate", fitness, X=X_train, Y=Y_train)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
    toolbox.register("select", tools.selTournament, tournsize=10)
    toolbox.register("gen_idx", gen, X_train.shape[1])
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.gen_idx, output_feature_number)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual, population)
    population = toolbox.population()
    halloffame = tools.HallOfFame(1)
    
    pop, logbook = algorithms.eaSimple(population, toolbox, cxpb, mutpb, ngen, halloffame=halloffame, verbose=False)
    
    fmask = halloffame[0]
    filtered_X_train = X_train[:, fmask]
    filtered_X_test = X_test[:, fmask]
    
    filtered_fnames = fnames[fmask]
    train_result = validate(X=filtered_X_train, Y=Y_train, task="train", model=model, **model_p)
    loo_result = validate(X=filtered_X_train, Y=Y_train, task="loo", model=model, **model_p)
    test_result = validate(X=filtered_X_train, Y=Y_train, task="test", model=model, **model_p, Xtest=filtered_X_test, Ytest=Y_test)
    
    print(f"best pop: {halloffame[0]}\n")
    print(f"train R2 is {train_result[1]} train RMSE is {train_result[0]}")
    print(f"loo R2 is {loo_result[1]} loo RMSE is {loo_result[0]}")
    print(f"test R2 is {test_result[1]} test RMSE is {test_result[0]}")
    
    filter_train = train.to_df().iloc[:, [True]+fmask]
    filter_test = test.to_df().iloc[:, [True]+fmask]
    # filter_test.to_csv("filtered_test.csv")
    # filter_train.to_csv("filtered_train.csv")
    
    from pathlib import Path
    logs = Path("ga-gbr-logs2.csv")
    
    if not logs.exists():
        with open(logs, "a") as f:
            f.writelines("random_state,individuals,length,loo_rmse,loo_r2,test_rmse,test_r2\n")
    if loo_result[1] > 0.6 and test_result[1] > 0.6 and loo_result[0] < test_result[0]:
        with open(logs, "a") as f:
            f.writelines(f"{random_state},{'-'.join([str(i) for i in fmask])},{len(set(fmask))},{loo_result[0]},{loo_result[1]},{test_result[0]},{test_result[1]}\n")

