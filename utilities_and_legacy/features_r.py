# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 12:24:54 2023

@author: senji
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 读取训练数据文件
train_data = pd.read_csv("train.csv")

# 提取特征部分
features = train_data.iloc[:, 2:]

# 计算皮尔逊相关系数矩阵
correlation_matrix = features.corr()

# 设置绘图风格
sns.set(style="white")

# 绘制热图并保存为高清图片
plt.figure(figsize=(10, 8))
sns.heatmap(correlation_matrix, annot=False, cmap="coolwarm", center=0)
plt.title("Pearson Correlation Heatmap")
plt.savefig("correlation_heatmap.png", dpi=300)  # 保存为高清图片
plt.show()

# 打印相关系数矩阵
print("Pearson Correlation Coefficients:")
print(correlation_matrix)

# 将相关系数矩阵保存为CSV文件
correlation_matrix.to_csv("correlation_matrix.csv", index=False)
