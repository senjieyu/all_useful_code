

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 11:25:10 2023

@author: senji
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_predict, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error, mean_absolute_error,r2_score
import pickle

# 读取CSV文件
data = pd.read_csv('8.18mrmr前20，RF前7_train.csv')

# 提取特征和目标变量
X = data.iloc[:, 2:].values
y = data.iloc[:, 1].values

# 特征标准化
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 创建SVR模型
svr = SVR(kernel='rbf')

# 设置交叉验证
kf = KFold(n_splits=5, shuffle=True, random_state=42)

# 初始化指标列表
pearson_scores = []  # 用于存储皮尔逊相关系数
rmse_scores = []
mae_scores = []
r_scores = []
# 进行交叉验证
for train_index, test_index in kf.split(X_scaled):
    X_train, X_test = X_scaled[train_index], X_scaled[test_index]
    y_train, y_test = y[train_index], y[test_index]

    # 训练模型
    svr.fit(X_train, y_train)

    # 预测
    y_pred = svr.predict(X_test)

    # 计算皮尔逊相关系数
    pearson_scores.append(np.corrcoef(y_test, y_pred)[0, 1])

    # 计算指标
    rmse_scores.append(np.sqrt(mean_squared_error(y_test, y_pred)))
    mae_scores.append(mean_absolute_error(y_test, y_pred))
    r_scores.append(r2_score(y_test, y_pred))
# 输出交叉验证指标
print("Pearson Correlation Coefficients:", pearson_scores)
print("RMSE Scores:", rmse_scores)
print("MAE Scores:", mae_scores)
print("R2 Scores:", r_scores)
# 计算平均指标
print("Average Pearson Correlation Coefficient:", np.mean(pearson_scores))
print("Average RMSE Score:", np.mean(rmse_scores))
print("Average MAE Score:", np.mean(mae_scores))
print("Average R2 Score:", np.mean(r_scores))

kernel_used = svr.kernel
print("Kernel used:", kernel_used)
gamma_used = svr.gamma
print("Gamma used:", gamma_used)
C_used = svr.C
print("C used:", C_used)


plt.figure(figsize=(10, 6))
plt.scatter(y_test, y_pred, c='blue', alpha=0.7)  # 使用y_test代替y
plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], linestyle='--', color='red', linewidth=2)
plt.title('Actual vs. Predicted')
plt.xlabel('Actual Values')
plt.ylabel('Predicted Values')
plt.show()
# 保存训练好的模型为pkl文件
with open('svr_model_8.18前7_features.pkl', 'wb') as f:
    pickle.dump(svr, f)



