# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 12:47:20 2023

@author: senji
"""

'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 读取训练数据文件
train_data = pd.read_csv("train.csv")

# 提取特征部分
features = train_data.iloc[:, 2:]

# 计算皮尔逊相关系数矩阵
correlation_matrix = features.corr()

# 找出相关系数大于0.95的特征对
high_correlation_pairs = (correlation_matrix.abs() > 0.95) & (correlation_matrix != 1)

# 筛选出相关系数大于0.95的特征名
high_correlation_features = set()
for feature in high_correlation_pairs.columns:
    correlated_features = list(high_correlation_pairs[feature][high_correlation_pairs[feature]].index)
    high_correlation_features.update(correlated_features)

# 提取相关系数大于0.95的特征对的相关系数矩阵
high_correlation_matrix = correlation_matrix.loc[high_correlation_features, high_correlation_features]

# 设置绘图风格
sns.set(style="white")

# 绘制相关系数大于0.95的特征对的热图并保存为高清图片
plt.figure(figsize=(10, 8))
sns.heatmap(high_correlation_matrix, annot=False, cmap="coolwarm", center=0)
plt.title("High Correlation Features Heatmap")
#plt.savefig("high_correlation_heatmap.png", dpi=300)  # 保存为高清图片
plt.show()
'''


#打印相关系数高于0.95的特征对
import pandas as pd

# 读取训练数据文件
train_data = pd.read_csv("train.csv")

# 提取特征部分
features = train_data.iloc[:, 2:]

# 计算皮尔逊相关系数矩阵
correlation_matrix = features.corr()

# 找出相关系数大于0.95的特征对
high_correlation_pairs = (correlation_matrix.abs() > 0.9) & (correlation_matrix != 1)

# 用集合来存储已经处理过的特征对，避免重复打印
processed_pairs = set()

# 打印相关系数大于0.95的特征对
for feature in high_correlation_pairs.columns:
    correlated_features = list(high_correlation_pairs[feature][high_correlation_pairs[feature]].index)
    for correlated_feature in correlated_features:
        # 判断是否已经处理过这个特征对，或者其对称特征对是否已经处理过
        if (correlated_feature, feature) not in processed_pairs:
            correlation_value = correlation_matrix.loc[feature, correlated_feature]
            print(f"{feature} : {correlated_feature} = {correlation_value:.2f}")
            processed_pairs.add((feature, correlated_feature))
            

'''
#大于0.9的热图
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 读取训练数据文件
train_data = pd.read_csv("train.csv")

# 提取特征部分
features = train_data.iloc[:, 2:]

# 计算皮尔逊相关系数矩阵
correlation_matrix = features.corr()

# 找出相关系数大于0.95的特征对
high_correlation_pairs = (correlation_matrix.abs() > 0.9) & (correlation_matrix != 1)

# 创建一个新的相关系数矩阵，只包含相关系数大于0.95的特征对
high_correlation_matrix = pd.DataFrame(index=features.columns, columns=features.columns, dtype=float)

# 用集合来存储已经处理过的特征对，避免重复填充
processed_pairs = set()

# 填充相关系数矩阵
for feature in high_correlation_pairs.columns:
    correlated_features = list(high_correlation_pairs[feature][high_correlation_pairs[feature]].index)
    for correlated_feature in correlated_features:
        # 判断是否已经处理过这个特征对，或者其对称特征对是否已经处理过
        if (correlated_feature, feature) not in processed_pairs:
            correlation_value = correlation_matrix.loc[feature, correlated_feature]
            high_correlation_matrix.loc[feature, correlated_feature] = correlation_value
            processed_pairs.add((feature, correlated_feature))

# 删除所有缺失值（无相关系数的位置）
high_correlation_matrix.dropna(axis=0, how='all', inplace=True)
high_correlation_matrix.dropna(axis=1, how='all', inplace=True)

# 设置绘图风格
sns.set(style="white")

# 绘制相关系数大于0.95的特征对的热图并保存为高清图片
plt.figure(figsize=(10, 8))
sns.heatmap(high_correlation_matrix, annot=True, cmap="coolwarm", center=0)
plt.title("High Correlation Features Heatmap")
plt.savefig("high_correlation_heatmap.png", dpi=300)  # 保存为高清图片
plt.show()   '''



