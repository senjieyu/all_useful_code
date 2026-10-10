# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 11:27:33 2024

@author: senji
"""

'''joblib 是一个用于在 Python 中进行高效的对象序列化和反序列化的库。
它通常用于在大规模数据处理和机器学习任务中保存和加载模型、数据等。以下是一些基本用法示例：




基本用法：
保存对象到文件：
import joblib
# 示例对象
data = {'name': 'John', 'age': 30, 'city': 'New York'}
# 保存到文件
joblib.dump(data, 'example.joblib')




从文件加载对象：
import joblib
# 从文件加载对象
loaded_data = joblib.load('example.joblib')
# 打印加载后的对象
print(loaded_data)



在机器学习中的应用：
from sklearn.datasets import load_iris
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import joblib
# 加载鸢尾花数据集
iris = load_iris()
X, y = iris.data, iris.target
# 划分数据集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
# 构建随机森林分类器
clf = RandomForestClassifier(n_estimators=100, random_state=42)
# 训练模型
clf.fit(X_train, y_train)
# 保存模型到文件
joblib.dump(clf, 'random_forest_model.joblib')
# 加载模型
loaded_model = joblib.load('random_forest_model.joblib')
# 在测试集上进行预测
predictions = loaded_model.predict(X_test)
# 打印预测结果
print(predictions)
这个例子中，我们使用 joblib 来保存和加载随机森林模型。在实际应用中，你可以替换模型和数据为你自己的模型和数据。joblib 的一个优点是它对于处理大型 NumPy 数组和复杂的数据结构时非常高效。'''