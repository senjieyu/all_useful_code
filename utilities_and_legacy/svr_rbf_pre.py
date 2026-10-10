# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 11:46:33 2023

@author: senji
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import pickle

# 读取CSV文件
data = pd.read_csv('selected_7_features_test.csv')

# 提取特征和目标变量
X = data.iloc[:, 2:].values
y = data.iloc[:, 1].values

# 特征标准化（使用之前训练集上的标准化器）
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 加载保存的SVR模型
with open('svr_model_7_features.pkl', 'rb') as f:
    svr = pickle.load(f)

# 预测测试集
y_pred = svr.predict(X_scaled)

# 计算评价指标
rmse = np.sqrt(mean_squared_error(y, y_pred))
mae = mean_absolute_error(y, y_pred)
r2 = r2_score(y, y_pred)

# 计算皮尔逊相关系数
pearson = np.corrcoef(y, y_pred)[0, 1]

# 输出评价指标
print("RMSE Score on Test Set:", rmse)
print("MAE Score on Test Set:", mae)
print("R2 Score on Test Set:", r2)
print("Pearson Correlation Coefficient on Test Set:", pearson)
