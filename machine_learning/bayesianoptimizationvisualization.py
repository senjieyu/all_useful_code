# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 17:47:17 2024

@author: senji
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

# 读取CSV文件
df = pd.read_excel('C:/Users/senji/Desktop/工作簿1.xlsx')  # 替换为你的文件路径

# 分离特征和目标变量
y = df.iloc[:, 0]  # 假设第一列是目标变量Bandgap
X = df.iloc[:, 1:]  # 其余列作为特征

# 应用PCA，计算第一主成分
pca = PCA(n_components=1)
X_pca = pca.fit_transform(X)

# 绘制散点图，其中X_pca是PCA1值
plt.figure(figsize=(10, 6))
scatter = plt.scatter(X_pca[:, 0], y, c=y, cmap='viridis', vmin=y.min(), vmax=y.max())
plt.colorbar(scatter, label='Bandgap (eV)')

# 添加拟合线
sns.regplot(x=X_pca[:, 0], y=y, scatter=False, color='orange')

plt.title('Bandgap vs. PCA1 After a Few Iterations')
plt.xlabel('PCA 1')
plt.ylabel('Bandgap (eV)')
plt.grid(True)
plt.show()
