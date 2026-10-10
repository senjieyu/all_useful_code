# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 11:21:02 2024

@author: senji
"""

# 导入所需的库和模块
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.datasets import load_iris
from sklearn.metrics import accuracy_score

# 加载鸢尾花数据集
iris = load_iris()
X, y = iris.data, iris.target

# 划分数据集为训练集和测试集（9:1的比例）
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)

# 构建机器学习流程
pipeline = Pipeline([
    ('scaler', StandardScaler()),  # 步骤1: 0-1标准化
    ('classifier', LogisticRegression())  # 步骤2: 逻辑回归
])

# 训练模型
pipeline.fit(X_train, y_train)

# 在测试集上进行预测
y_pred = pipeline.predict(X_test)

# 计算准确性
accuracy = accuracy_score(y_test, y_pred)

print(pipeline.score(X_train,y_train))
print(f'模型准确性: {accuracy:.2f}')
