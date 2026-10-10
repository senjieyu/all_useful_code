# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 15:31:54 2025

@author: senji
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jan 5 14:47:21 2025
Updated by Assistant
"""

import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_predict, KFold
from sklearn.cross_decomposition import PLSRegression
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.metrics import mean_squared_error
import numpy as np
import pandas as pd

# 读取CSV文件
df = pd.read_csv(
    r'C:/Users/senji/Desktop/高光谱及高光谱仅部位/高光谱定量模型/数据集/原始平滑.csv', 
    encoding='latin1'
)

# 提取目标值列和特征列
X = df.iloc[:, 9:]  # 假设特征列从第9列开始
y = df.iloc[:, 7]   # 假设目标值列为第2列

# 初始化最佳结果变量
min_rmse = np.inf
best_regressor = None
best_regressor_instance = None
best_dimensions = 0
best_X_train = None
best_X_test = None
best_y_train = None
best_y_test = None
best_y_pred = None
best_train_rmse =  0
best_tese_rmse =  0
best_train_r =  0
best_tese_r =  0
# 保存每个回归器的结果
results = {}
regressors_results = {}

# 定义回归模型
regressors = {
    "SVR": lambda: SVR(),
    "Linear Regression": lambda: LinearRegression(),
    "Ridge": lambda: Ridge(),
    "Lasso": lambda: Lasso(),
    "MLP": lambda: MLPRegressor(max_iter=500),
    "Random Forest": lambda: RandomForestRegressor(),
    "Decision Tree": lambda: DecisionTreeRegressor(),
    "KNN": lambda: KNeighborsRegressor()
}

# 循环遍历所有回归模型
for regressor_name, regressor_func in regressors.items():
    regressors_results[regressor_name] = {
        "dimensions": [],
        "train_rmse": [],
        "test_rmse": [],
        "train_corr": [],
        "test_corr": []
    }
    print(f"Processing {regressor_name}...")

    regressor = regressor_func()  # 动态实例化

    for n_components in range(2,30):  # n_components 不超过特征数量
        pls = PLSRegression(n_components=n_components)

        # 划分训练集和测试集
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        # 检查并处理缺失值
        if X_train.isnull().any().any() or X_test.isnull().any().any():
            X_train = X_train.fillna(X_train.mean())
            X_test = X_test.fillna(X_test.mean())

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

        # 计算皮尔逊相关系数
        train_corr = np.corrcoef(y_train, y_train_pred)[0, 1]
        test_corr = np.corrcoef(y_test, y_test_pred)[0, 1]

        # 保存每一维度的结果
        regressors_results[regressor_name]["dimensions"].append(n_components)
        regressors_results[regressor_name]["train_rmse"].append(train_rmse)
        regressors_results[regressor_name]["test_rmse"].append(test_rmse)
        regressors_results[regressor_name]["train_corr"].append(train_corr)
        regressors_results[regressor_name]["test_corr"].append(test_corr)

        # 保存最佳结果
        if test_rmse < min_rmse:
            min_rmse = test_rmse
            best_regressor = regressor_name  # 保存最佳模型名字
            best_regressor_instance = regressor  # 保存最佳模型实例
            best_dimensions = n_components
            best_X_train = X_train_pls
            best_X_test = X_test_pls
            best_y_train = y_train
            best_y_test = y_test
            best_y_pred = y_test_pred
            best_train_rmse =  train_rmse
            best_tese_rmse =  test_rmse
            best_train_r =  train_corr
            best_tese_r =  test_corr


        # 输出结果
        print(f"Regressor: {regressor_name}, Dimensions: {n_components}, Train RMSE: {train_rmse}, Test RMSE: {test_rmse}, Train Pearson Correlation: {train_corr}, Test Pearson Correlation: {test_corr}")

print('========================================================================')
print(f"\nBest Dimensionality Reduction Result:")
print(f"Regressor: {best_regressor}")
print(f"Dimensions: {best_dimensions}")
print(f"\nBest Dimensionality Reduction Result:")
print(f"Regressor: {best_regressor}")
print(f"Dimensions: {best_dimensions}")
print(f"Dimensions: {best_dimensions}")
print(f"best_train_rmse:{best_train_rmse}")
print(f"best_tese_rmse:{best_tese_rmse}")
print(f"best_train_r:{best_train_r}")
print(f"best_tese_r:{best_tese_r}")
# 保存最佳模型数据
best_X_train = pd.DataFrame(best_X_train).reset_index(drop=True)
best_y_train = pd.Series(best_y_train).reset_index(drop=True)
best_X_test = pd.DataFrame(best_X_test).reset_index(drop=True)
best_y_test = pd.Series(best_y_test).reset_index(drop=True)


train_data = pd.DataFrame(best_X_train)
train_data['Target'] = best_y_train
train_data.to_csv(r'C:/Users/senji/Desktop/高光谱及高光谱仅部位/高光谱定量模型/代码/最佳模型数据集//best_train_data.csv', index=False, encoding='utf-8')

test_data = pd.DataFrame(best_X_test)
test_data['Target'] = best_y_test
test_data.to_csv(r'C:/Users/senji/Desktop/高光谱及高光谱仅部位/高光谱定量模型/代码/最佳模型数据集//best_test_data.csv', index=False, encoding='utf-8')

# 从文件加载训练集和测试集
train_data = pd.read_csv(r'C:/Users/senji/Desktop/高光谱及高光谱仅部位/高光谱定量模型/代码/最佳模型数据集//best_train_data.csv')
test_data = pd.read_csv(r'C:/Users/senji/Desktop/高光谱及高光谱仅部位/高光谱定量模型/代码/最佳模型数据集//best_test_data.csv')

X_train_loaded = train_data.iloc[:, :-1]  # 特征列
y_train_loaded = train_data['Target']    # 目标值列
X_test_loaded = test_data.iloc[:, :-1]   # 特征列
y_test_loaded = test_data['Target']      # 目标值列

# 动态重新实例化最佳模型
best_regressor_instance = regressors[best_regressor]()  # 动态实例化

# 重新训练模型
best_regressor_instance.fit(X_train_loaded, y_train_loaded)

# 预测
y_train_pred_loaded = best_regressor_instance.predict(X_train_loaded)
y_test_pred_loaded = best_regressor_instance.predict(X_test_loaded)

# 计算误差和交叉验证
train_relative_error_loaded = np.mean(np.abs((y_train_loaded - y_train_pred_loaded) / (y_train_loaded + 1e-10))) * 100
test_relative_error_loaded = np.mean(np.abs((y_test_loaded - y_test_pred_loaded) / (y_test_loaded + 1e-10))) * 100

kf_5 = KFold(n_splits=5, shuffle=True, random_state=42)
y_pred_5_fold = cross_val_predict(best_regressor_instance, X_train_loaded, y_train_loaded, cv=kf_5)
cross_val_5_relative_error_loaded = np.mean(np.abs((y_train_loaded - y_pred_5_fold) / (y_train_loaded + 1e-10))) * 100

kf_10 = KFold(n_splits=10, shuffle=True, random_state=42)
y_pred_10_fold = cross_val_predict(best_regressor_instance, X_train_loaded, y_train_loaded, cv=kf_10)
cross_val_10_relative_error_loaded = np.mean(np.abs((y_train_loaded - y_pred_10_fold) / (y_train_loaded + 1e-10))) * 100

# 打印误差结果
print("\nUsing Saved Train and Test Data:")
print(f"Train Relative Error (Loaded): {train_relative_error_loaded:.2f}%")
print(f"Test Relative Error (Loaded): {test_relative_error_loaded:.2f}%")
print(f"5-Fold Cross-Validation Relative Error (Loaded): {cross_val_5_relative_error_loaded:.2f}%")
print(f"10-Fold Cross-Validation Relative Error (Loaded): {cross_val_10_relative_error_loaded:.2f}%")

# 可视化RMSE和相关系数
plt.rcParams.update({'font.size': 12, 'figure.dpi': 150})
fig, axes = plt.subplots(len(regressors), 2, figsize=(18, len(regressors) * 6))

colors = ['#8ECFC9', '#FFBE7A', '#FA7F6F', '#82B0D2']

for idx, (regressor_name, result) in enumerate(regressors_results.items()):
    # 绘制训练集和测试集RMSE
    axes[idx, 0].plot(result["dimensions"], result["train_rmse"], marker='o', color=colors[0], label='Train RMSE')
    axes[idx, 0].plot(result["dimensions"], result["test_rmse"], marker='o', color=colors[1], label='Test RMSE')
    axes[idx, 0].set_xlabel('Dimensions', fontsize=14)
    axes[idx, 0].set_ylabel('RMSE', fontsize=14)
    axes[idx, 0].set_title(f'{regressor_name} RMSE vs Dimensions', fontsize=16)
    axes[idx, 0].legend(fontsize=12)

    # 绘制相关系数
    axes[idx, 1].plot(result["dimensions"], result["train_corr"], marker='o', color=colors[2], label='Train Pearson Correlation')
    axes[idx, 1].plot(result["dimensions"], result["test_corr"], marker='o', color=colors[3], label='Test Pearson Correlation')
    axes[idx, 1].set_xlabel('Dimensions', fontsize=14)
    axes[idx, 1].set_ylabel('Pearson Correlation', fontsize=14)
    axes[idx, 1].set_title(f'{regressor_name} Correlation vs Dimensions', fontsize=16)
    axes[idx, 1].legend(fontsize=12)

plt.tight_layout()
plt.show()