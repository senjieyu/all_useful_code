# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 16:53:05 2024

@author: senji
"""

# =============================================================================
# # 连续投影降维
# 
# import pandas as pd
# from sklearn.model_selection import train_test_split
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.metrics import accuracy_score, confusion_matrix
# 
# # 读取训练集和测试集
# df_train = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林训练集.csv')
# df_test = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林测试集.csv')
# 
# # 提取目标值列和特征列
# X_train = df_train.iloc[:, 1:]  # 假设特征列从第三列开始
# y_train = df_train.iloc[:, 0]   # 假设目标值列为第二列
# 
# X_test = df_test.iloc[:, 1:]  # 假设特征列从第三列开始
# y_test = df_test.iloc[:, 0]   # 假设目标值列为第二列
# 
# 
# 
# 
# 
# # 初始化随机森林分类器
# random_forest_classifier = RandomForestClassifier()
# 
# # 训练分类器模型
# random_forest_classifier.fit(X_train, y_train)
# 
# # 预测训练集和测试集
# train_predictions = random_forest_classifier.predict(X_train)
# test_predictions = random_forest_classifier.predict(X_test)
# 
# # 计算训练集和测试集准确率
# train_accuracy = accuracy_score(y_train, train_predictions)
# test_accuracy = accuracy_score(y_test, test_predictions)
# 
# 
# # 输出结果
# print(f"Train Accuracy: {train_accuracy}, Test Accuracy: {test_accuracy}")
# 
# =============================================================================





# =============================================================================
# #随机网格搜索
# import pandas as pd
# from sklearn.model_selection import train_test_split, RandomizedSearchCV
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.metrics import accuracy_score
# from sklearn.model_selection import train_test_split
# from sklearn.model_selection import cross_val_score, LeaveOneOut
# from sklearn.ensemble import RandomForestClassifier
# 
# # 读取训练集和测试集
# # 读取训练集和测试集
# df_train = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林训练集.csv')
# df_test = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林测试集.csv')
# 
# # 合并两个数据集
# df_combined = pd.concat([df_train, df_test], ignore_index=True)
# X = df_combined.iloc[:, 1:]  # 假设特征列从第二列开始
# y = df_combined.iloc[:, 0]   # 假设目标值列为第一列
# 
# 
# 
# # 将数据集随机划分成训练集和测试集，按照4:1的比例
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
# 
# 
# 
# # 提取特征列和目标值列
# 
# 
# 
# # 定义随机森林分类器
# random_forest_classifier = RandomForestClassifier()
# 
# # 定义超参数搜索范围
# param_distributions = {
#     'n_estimators': [50, 100, 150, 200, 250, 300],
#     'max_depth': [None, 5, 10, 15, 20, 25],
#     'min_samples_split': [2, 3, 4, 5, 6, 7, 8, 9, 10],
#     'min_samples_leaf': [1, 2, 3, 4, 5,6]
# }
# 
# # 使用RandomizedSearchCV进行随机搜索
# random_search = RandomizedSearchCV(estimator=random_forest_classifier, param_distributions=param_distributions, n_iter=10, cv=5)
# random_search.fit(X_train, y_train)
# 
# # 获取最佳参数
# best_params = random_search.best_params_
# 
# # 使用最佳参数的分类器重新训练模型
# best_classifier = RandomForestClassifier(**best_params)
# 
# # 使用LOOCV评估模型
# loo = LeaveOneOut()
# accuracies = cross_val_score(best_classifier, X_train, y_train, cv=loo)
# 
# # 计算 LOOCV 的平均准确率
# average_accuracy = accuracies.mean()
# 
# # 输出结果
# print(f"Best Parameters: {best_params}")
# print(f"LOOCV Accuracy: {average_accuracy}")
# =============================================================================

# =============================================================================
# # # 网格搜索
# # import pandas as pd
# # from sklearn.model_selection import train_test_split, GridSearchCV
# # from sklearn.ensemble import RandomForestClassifier
# # from sklearn.metrics import accuracy_score
# # from sklearn.model_selection import train_test_split
# # from sklearn.model_selection import cross_val_score, LeaveOneOut
# # from sklearn.ensemble import RandomForestClassifier
# 
# # train_data = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林训练集.csv')
# # test_data = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林测试集.csv')
# 
# # # 提取特征列和目标列
# # X_train = train_data.iloc[:, 1:]  # 特征列，从第二列开始到最后一列
# # y_train = train_data.iloc[:, 0]   # 目标列，第一列
# 
# # X_test = test_data.iloc[:, 1:]
# # y_test = test_data.iloc[:, 0]
# 
# # # 定义随机森林分类器
# # random_forest_classifier = RandomForestClassifier(random_state=42,n_jobs=12)
# 
# 
# 
# # # 定义超参数搜索范围
# # param_grid = {
# #     'n_estimators': [5,10, 15,30, 50, 100, 150, 200, 250, ],
# #     'max_depth': [None, 5, 10, 15, 20, 25],
# #     'min_samples_split': [ 2, 3, 4, 5, 6, 7, 8, 9, 10],
# #     'min_samples_leaf': [1, 2, 3, 4, 5, 6], 
# #     "max_features": ["log2","sqrt",None]
# # }
# 
# # # 使用GridSearchCV进行网格搜索
# # grid_search = GridSearchCV(estimator=random_forest_classifier, param_grid=param_grid, cv=5)
# # grid_search.fit(X_train, y_train)
# 
# # # 获取最佳参数
# # best_params = grid_search.best_params_
# 
# # # 使用最佳参数的分类器重新训练模型
# # best_classifier = RandomForestClassifier(**best_params)
# 
# # # 使用LOOCV评估模型
# # loo = LeaveOneOut()
# # accuracies = cross_val_score(best_classifier, X_train, y_train, cv=5)
# 
# # # 计算 LOOCV 的平均准确率
# # average_accuracy = accuracies.mean()
# 
# # # 输出结果
# # print(f"Best Parameters: {best_params}")
# # print(f"LOOCV Accuracy: {average_accuracy}")
# =============================================================================




# loocv 网格搜索
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score, LeaveOneOut
from sklearn.ensemble import RandomForestClassifier

train_data = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林训练集.csv')
test_data = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林测试集.csv')

# 提取特征列和目标列
X_train = train_data.iloc[:, 1:]  # 特征列，从第二列开始到最后一列
y_train = train_data.iloc[:, 0]   # 目标列，第一列

X_test = test_data.iloc[:, 1:]
y_test = test_data.iloc[:, 0]

# 定义随机森林分类器
random_forest_classifier = RandomForestClassifier(random_state=42,n_jobs=24)

# 定义超参数搜索范围
param_grid = {
    'n_estimators': [100,150,200],
    'max_depth': [None, 2,4,6],
    'min_samples_split': [2, 3, 4, 5],
    # 'min_samples_leaf': [1, 2, 3],
    'random_state' : list(range(1, 101, 1))
}

# 使用GridSearchCV进行网格搜索
grid_search = GridSearchCV(estimator=random_forest_classifier, param_grid=param_grid, cv=5)
grid_search.fit(X_train, y_train)

# 获取最佳参数
best_params = grid_search.best_params_

# 使用最佳参数的分类器重新训练模型
best_classifier = RandomForestClassifier(**best_params)

# 使用LOOCV评估模型
loo = LeaveOneOut()
accuracies = cross_val_score(best_classifier, X_train, y_train, cv=5)

# 计算 LOOCV 的平均准确率
average_accuracy = accuracies.mean()

# 输出结果
print(f"Best Parameters: {best_params}")
print(f"5CV Accuracy: {average_accuracy}")


# =============================================================================
# Best Parameters: {'max_depth': 5, 'min_samples_leaf': 1, 'min_samples_split': 3, 'n_estimators': 100}
# LOOCV Accuracy: 0.8846153846153846
# =============================================================================




