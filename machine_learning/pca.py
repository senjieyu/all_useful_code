# -*- coding: utf-8 -*-
"""
Created on Thu May 30 10:16:44 2024

@author: senji
"""

import pandas as pd, numpy as np, re
import matplotlib.pyplot as plt
plt.rcParams["figure.dpi"] = 300
plt.rcParams["font.sans-serif"] = ["SimHei"]
plt.rcParams["font.size"] = 12
plt.rcParams["axes.unicode_minus"] = False
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.colors import Normalize

data = pd.read_excel(r'C:/Users/senji/Desktop/低熔点合金突破文章/全部数据集--用于模式识别.xlsx' )
data.index = data.index.astype(str)

X = data.iloc[:, 2:]
Y = data.iloc[:, 1].to_frame()


plot_X = StandardScaler().fit_transform(X)
plot_X = PCA().fit_transform(plot_X)
fig, ax = plt.subplots()
norm = Normalize(vmin=Y.min(), vmax=Y.max())
scatter = ax.scatter(plot_X[:, 0], plot_X[:, 1], c=Y.values, cmap="coolwarm", norm=norm)
cb = fig.colorbar(scatter, extend='both')
cb.ax.set_title("COP")
ax.set_title("删除融冰/蓄冰后的数据分布")
ax.set_xlabel("PCA1")
ax.set_ylabel("PCA2")




# 筛选出Y小于30的点
mask = Y.values.flatten() < 30
X_filtered = X[mask]
Y_filtered = Y[mask]

# 对筛选后的数据进行标准化和PCA降维
plot_X_filtered = StandardScaler().fit_transform(X_filtered)
plot_X_filtered = PCA().fit_transform(plot_X_filtered)

# 绘制筛选后的数据点
fig_filtered, ax_filtered = plt.subplots()
norm_filtered = Normalize(vmin=Y_filtered.min(), vmax=Y_filtered.max())
scatter_filtered = ax_filtered.scatter(plot_X_filtered[:, 0], plot_X_filtered[:, 1], c=Y_filtered.values, cmap="coolwarm", norm=norm_filtered)
cb_filtered = fig_filtered.colorbar(scatter_filtered, extend='both')
cb_filtered.ax.set_title("COP")
ax_filtered.set_title("Y小于30的点的数据分布")
ax_filtered.set_xlabel("PCA1")
ax_filtered.set_ylabel("PCA2")

plt.show()