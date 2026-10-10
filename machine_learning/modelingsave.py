# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:17:54 2024

@author: senji
"""

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.cross_decomposition import PLSRegression
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import accuracy_score

df = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/部位一阶导函数.csv')

# 提取目标值列和特征列
X = df.iloc[:, 2:]  # 假设特征列从第三列开始
y = df.iloc[:, 1]   # 假设目标值列为第二列

# 初始化最佳结果变量
max_accuracy = 0
best_classifier = None
best_dimensions = 0
best_X_train_pls = None
best_X_test_pls = None
best_y_train = None
best_y_test = None
best_y_pred = None

classifiers = {
    "PLS_SVM": SVC(),
    "Logistic Regression": LogisticRegression(max_iter=1000,random_state=42),
    "K-Nearest Neighbors": KNeighborsClassifier(),
    "Decision Tree": DecisionTreeClassifier(random_state=42),
    "Random Forest": RandomForestClassifier(random_state=42),
    "Gradient Boosting": GradientBoostingClassifier(random_state=42)
}

for classifier_name, classifier in classifiers.items():
    for n_components in range(10, 25):  # 尝试2到100个维度
        pls = PLSRegression(n_components=n_components)

        # 划分训练集和测试集
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        # 训练PLS模型
        pls.fit(X_train, y_train)

        # 使用训练好的模型对数据进行降维
        X_train_pls = pls.transform(X_train)
        X_test_pls = pls.transform(X_test)

        # 训练分类器模型
        classifier.fit(X_train_pls, y_train)

        # 预测测试集
        y_pred = classifier.predict(X_test_pls)
        train_accuracy = accuracy_score(y_train, classifier.predict(X_train_pls))
        # 计算测试集准确率
        test_accuracy = accuracy_score(y_test, y_pred)
        print(f"Classifier: {classifier_name}, Dimensions: {n_components}, Train Accuracy: {train_accuracy}, Test Accuracy: {test_accuracy}")
        # 保存最佳结果数据集、训练集和测试集
        if test_accuracy > max_accuracy:
            max_accuracy = test_accuracy
            best_classifier = classifier_name
            best_dimensions = n_components
            best_X_train_pls = X_train_pls
            best_X_test_pls = X_test_pls
            best_y_train = y_train
            best_y_test = y_test
            best_y_pred = y_pred

# 打印最佳降维结果
print(f"Best Dimensionality Reduction Result - Dimensions: {best_dimensions}, Classifier: {best_classifier}")

# 保存最佳结果数据集、训练集和测试集
# 保存最佳结果数据集、训练集和测试集
best_X_train_pls_df = pd.DataFrame(best_X_train_pls, columns=[f"PLS_Component_{i+1}" for i in range(best_dimensions)])
best_X_test_pls_df = pd.DataFrame(best_X_test_pls, columns=[f"PLS_Component_{i+1}" for i in range(best_dimensions)])

# 合并特征数据和目标数据
combined_train_pls = pd.concat([best_y_train.reset_index(drop=True), best_X_train_pls_df], axis=1)
combined_test_pls = pd.concat([best_y_test.reset_index(drop=True), best_X_test_pls_df], axis=1)

# 保存合并后的训练集和测试集
combined_train_pls.to_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林训练集.csv', index=False)
combined_test_pls.to_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林测试集.csv', index=False)

