# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 20:34:02 2024

@author: senji
"""


from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import load_iris

# 加载鸢尾花数据集
iris = load_iris()
X, y = iris.data, iris.target

# 划分训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 使用StandardScaler进行Z-归一化
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# 打印归一化后的训练集和测试集
print("归一化后的训练集:")
print(X_train_scaled)
print("\n归一化后的测试集:")
print(X_test_scaled)
