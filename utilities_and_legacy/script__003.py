# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:12:51 2024

@author: senji
"""

# 导入必要的库
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import mean_squared_error
from sklearn.datasets import load_iris

# 加载鸢尾花数据集
iris = load_iris()
X = iris.data
y = iris.target

# 划分训练集和测试集（4:1）
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 数据归一化（Z-归一化）
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# 定义逻辑回归模型
logistic_regression = LogisticRegression(max_iter=int(1e6), solver='saga')

# 定义超参数网格
param_grid = {'C': [0.001, 0.01, 0.1, 1, 10, 100], 'penalty': ['l1', 'l2']}

# 使用网格搜索进行超参数优化，评估指标为RMSE，采用10折交叉验证
grid_search = GridSearchCV(logistic_regression, param_grid, cv=10, scoring='neg_mean_squared_error')
grid_search.fit(X_train_scaled, y_train)

# 打印最佳参数
best_params = grid_search.best_params_
print("最佳参数：", best_params)

# 在测试集上进行预测
y_pred = grid_search.predict(X_test_scaled)

# 计算RMSE
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
print("测试集上的RMSE：", rmse)
