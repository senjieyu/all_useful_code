# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:15:49 2024

@author: senji
"""

import numpy as np
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.datasets import load_iris
from sklearn.pipeline import Pipeline
from sklearn.multiclass import OneVsRestClassifier

# 加载鸢尾花数据集
iris = load_iris()
X = iris.data
y = iris.target

# 划分训练集和测试集（4:1）
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 创建Pipeline
pipeline = Pipeline([
    ('scaler', StandardScaler()),  # Z-归一化
    ('classifier', OneVsRestClassifier(LogisticRegression(max_iter=int(1e6), solver='saga')))  # 使用OneVsRest包装逻辑回归
])

# 定义超参数网格
param_grid = {'classifier__estimator__C': [0.001, 0.01, 0.1, 1, 10, 100], 'classifier__estimator__penalty': ['l1', 'l2']}

# 使用GridSearchCV进行超参数优化，评估指标为ROC-AUC，采用10折交叉验证
grid_search = GridSearchCV(pipeline, param_grid, cv=10, scoring='roc_auc_ovr', n_jobs=-1)  # 注意 scoring 参数   多分类的AUC-ROC
grid_search.fit(X_train, y_train)

# 在测试集上进行预测
y_pred_proba = grid_search.predict_proba(X_test)  # 预测类别概率
y_pred = np.argmax(y_pred_proba, axis=1)  # 获取最大概率对应的类别

# 计算ROC-AUC
roc_auc = roc_auc_score(y_test, y_pred_proba, multi_class='ovr')  # 使用One-vs-Rest的策略
print("测试集上的ROC-AUC：", roc_auc)

# 打印最佳参数
best_params = grid_search.best_params_
print("最佳参数：", best_params)

'''
AUC-ROC (Area Under the Receiver Operating Characteristic curve) 是一种用于评估二分类和多分类分类器性能的指标。虽然名字中包含 "Receiver Operating Characteristic"，但它不仅适用于二分类问题，也可以扩展到多分类问题。

AUC-ROC for Binary Classification:

在二分类问题中，AUC-ROC是一种度量分类器在不同阈值下真正例率（True Positive Rate，也称为灵敏度或召回率）和假正例率（False Positive Rate）之间的性能的方法。ROC曲线是以不同的分类阈值为基础，绘制出真正例率和假正例率之间的关系，AUC-ROC就是ROC曲线下的面积。

对于二分类问题，AUC-ROC的值介于0和1之间。完美的分类器AUC-ROC为1，随机分类器AUC-ROC为0.5。当AUC-ROC大于0.5时，说明模型优于随机猜测。

AUC-ROC for Multiclass Classification:

在多分类问题中，AUC-ROC的计算方式略有不同。通常，使用一种称为"one-vs-rest"（OvR）或"one-vs-all"的策略。对于每个类别，将其余所有类别视为一个类别，计算该类别与其他所有类别的AUC-ROC。然后，将所有类别的AUC-ROC取平均值作为多分类AUC-ROC。

在scikit-learn中，通常使用multi_class='ovr'参数来指定"one-vs-rest"的策略。因此，对于多分类问题，AUC-ROC的评估指标可以写作 roc_auc_ovr，表示"One-vs-Rest"的AUC-ROC。

总的来说，AUC-ROC是一种广泛用于评估分类器性能的指标，可以适用于二分类和多分类问题。'''