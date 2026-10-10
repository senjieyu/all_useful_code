# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 17:19:34 2024

@author: senji
"""

import pandas as pd
from sklearn.model_selection import train_test_split, LeaveOneOut, cross_val_predict
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC  # 以SVM分类器为例，你可以根据需求选择其他分类器
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

# 读取训练集和测试集csv文件
train_data = pd.read_csv(r'C:/Users/senji/Desktop/图表/最佳/111best_X_train_pls.csv')
test_data = pd.read_csv(r'C:/Users/senji/Desktop/图表/最佳/111best_X_test_pls.csv')

# 提取特征列和目标列
X_train = train_data.iloc[:, 1:]  # 特征列，从第二列开始到最后一列
y_train = train_data.iloc[:, 0]   # 目标列，第一列

X_test = test_data.iloc[:, 1:]
y_test = test_data.iloc[:, 0]

# 使用LabelEncoder将目标列转换为数值
label_encoder = LabelEncoder()
y_train = label_encoder.fit_transform(y_train)
y_test = label_encoder.transform(y_test)

# 初始化分类器
'''这里修改模型，只需要修改文件名和这里'''
# classifier = LogisticRegression(max_iter=1000,random_state=42)
# classifier = SVC(random_state=42)
# classifier = KNeighborsClassifier()
# classifier = DecisionTreeClassifier(random_state=42)
classifier = RandomForestClassifier(random_state=42)
# classifier = GradientBoostingClassifier(random_state=42)

# 拟合模型
classifier.fit(X_train, y_train)

# 预测训练集和测试集
train_predictions = classifier.predict(X_train)
test_predictions = classifier.predict(X_test)

# 计算混淆矩阵
train_confusion_matrix = confusion_matrix(y_train, train_predictions)
test_confusion_matrix = confusion_matrix(y_test, test_predictions)

print("训练集混淆矩阵:")
print(train_confusion_matrix)
print("\n测试集混淆矩阵:")
print(test_confusion_matrix)


