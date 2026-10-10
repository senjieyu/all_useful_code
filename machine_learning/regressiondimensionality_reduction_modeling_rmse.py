# -*- coding: utf-8 -*-
"""
Created on Tue May  7 16:05:01 2024

@author: senji
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.cross_decomposition import PLSRegression
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.metrics import mean_squared_error

# 读取CSV文件
df = pd.read_csv(r'C:/Users/senji/Desktop/梅菜扣肉/一阶导函数.csv')

# 提取目标值列和特征列
X = df.iloc[:, 2:]  # 假设特征列从第三列开始
y = df.iloc[:, 0]   # 假设目标值列为第二列

# 定义回归模型
regressors = {
    "SVR": SVR(),
    "Linear Regression": LinearRegression(),
    "Ridge": Ridge(),
    "Lasso": Lasso(),
    "MLP": MLPRegressor(),
    "Random Forest": RandomForestRegressor(),
    "Decision Tree": DecisionTreeRegressor(),
    "KNN": KNeighborsRegressor()
}

# Initialize lists to store RMSE values
train_rmse_values = {regressor_name: [] for regressor_name in regressors}
test_rmse_values = {regressor_name: [] for regressor_name in regressors}

# 初始化最佳结果变量
min_rmse = np.inf
best_regressor = None
best_dimensions = 0
best_X_train = None
best_X_test = None
best_y_train = None
best_y_test = None
best_y_pred = None

# 循环遍历所有回归模型
for regressor_name, regressor in regressors.items():
    for n_components in range(2, 30):  # 尝试2到30个维度
        pls = PLSRegression(n_components=n_components)

        # 划分训练集和测试集
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        # 训练PLS模型
        pls.fit(X_train, y_train)

        # 使用训练好的模型对数据进行降维
        X_train_pls = pls.transform(X_train)
        X_test_pls = pls.transform(X_test)

        # 训练回归模型
        regressor.fit(X_train_pls, y_train)

        # 预测训练集和测试集
        y_train_pred = regressor.predict(X_train_pls)
        y_test_pred = regressor.predict(X_test_pls)

        # 计算RMSE
        train_rmse = mean_squared_error(y_train, y_train_pred, squared=False)
        test_rmse = mean_squared_error(y_test, y_test_pred, squared=False)

        # 保存最佳结果
        if test_rmse < min_rmse:
            min_rmse = test_rmse
            best_regressor = regressor_name
            best_dimensions = n_components
            best_X_train = X_train_pls
            best_X_test = X_test_pls
            best_y_train = y_train
            best_y_test = y_test
            best_y_pred = y_test_pred

        # 保存RMSE到对应的列表中
        train_rmse_values[regressor_name].append(train_rmse)
        test_rmse_values[regressor_name].append(test_rmse)

        # 输出结果
        print(f"Regressor: {regressor_name}, Dimensions: {n_components}, Train RMSE: {train_rmse}, Test RMSE: {test_rmse}")

# 打印最佳降维结果
print(f"Best Dimensionality Reduction Result - Dimensions: {best_dimensions}, Regressor: {best_regressor}, Test RMSE: {min_rmse}")

# 绘制RMSE图像
# 绘制RMSE图像


from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.sans-serif'] = ['SimHei']  # 设置字体为黑体


plt.figure(figsize=(14, 6),dpi = 300)

# Plotting training set RMSE
plt.subplot(1, 2, 1)
for regressor_name, rmse_values in train_rmse_values.items():
    plt.plot(range(2, 30), rmse_values, label=regressor_name, marker='o')

plt.title('训练集 RMSE',fontproperties='SimHei', fontsize=16)
plt.xlabel('维度', fontproperties='SimHei',fontsize=16)
plt.ylabel('RMSE',fontproperties='SimHei', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend()
plt.grid(True)

# 绘制测试集RMSE
plt.subplot(1, 2, 2)
for regressor_name, rmse_values in test_rmse_values.items():
    plt.plot(range(2, 30), rmse_values, label=regressor_name, marker='o')

plt.title('测试集 RMSE', fontproperties='SimHei',fontsize=16)
plt.xlabel('维度',fontproperties='SimHei', fontsize=16)
plt.ylabel('RMSE',fontproperties='SimHei', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()


# 打印最佳降维结果
print(f"Best Dimensionality Reduction Result - Dimensions: {best_dimensions}, Regressor: {best_regressor}, Test RMSE: {min_rmse}")
