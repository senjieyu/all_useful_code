# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:13:08 2024

@author: senji
"""

import pandas as pd
from sklearn.model_selection import train_test_split, LeaveOneOut, cross_val_score
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC  # 以SVM分类器为例，你可以根据需求选择其他分类器
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import accuracy_score

# 读取训练集和测试集csv文件
train_data = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林训练集.csv')
test_data = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/随机森林测试集.csv')

# 提取特征列和目标列
X_train = train_data.iloc[:, 1:]  # 特征列，从第二列开始到最后一列
y_train = train_data.iloc[:, 0]   # 目标列，第一列

X_test = test_data.iloc[:, 1:]
y_test = test_data.iloc[:, 0]



# 初始化分类器
'''这里修改模型，只需要修改文件名和这里'''
# classifier = LogisticRegression ( max_iter=1000 )
# classifier = SVC()   
# classifier = KNeighborsClassifier()  
# classifier = DecisionTreeClassifier()  
classifier = RandomForestClassifier(random_state=42) 
# classifier = RandomForestClassifier(max_depth= 5, min_samples_leaf= 1, min_samples_split=3, n_estimators= 100,random_state=42)
# classifier = GradientBoostingClassifier(max_depth= 5, max_features='log2', min_samples_leaf=1, min_samples_split= 2, n_estimators= 250)




# 拟合模型并进行预测
classifier.fit(X_train, y_train)
train_predictions = classifier.predict(X_train)
test_predictions = classifier.predict(X_test)

# 计算并输出训练集和测试集的准确率
train_accuracy = accuracy_score(y_train, train_predictions)
test_accuracy = accuracy_score(y_test, test_predictions)

print("训练集准确率:", train_accuracy)
print("测试集准确率:", test_accuracy)

# 使用LOOCV计算并输出准确率
loo = LeaveOneOut()
cv_scores = cross_val_score(classifier, X_train, y_train, cv=loo)

print("LOOCV准确率:", cv_scores.mean())




