# -*- coding: utf-8 -*-
"""
Created on Fri May 31 11:13:00 2024

@author: senji
"""

# -*- coding: utf-8 -*-
"""
Created on Fri May 31 09:04:46 2024

@author: senji
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.colors import Normalize
from matplotlib import rcParams

# Set Times New Roman as the default font

rcParams['font.serif'] = ['Times New Roman']
plt.rcParams["figure.dpi"] = 300
plt.rcParams["font.sans-serif"] = ["Times New Roman"]
plt.rcParams["font.size"] = 12
plt.rcParams["axes.unicode_minus"] = False

# 读取数据
data = pd.read_csv(r'C:/Users/senji/Desktop/模式识别ANN/全部数据集--用于模式识别.csv')
data.index = data.index.astype(str)

# 提取特征和目标
X = data.iloc[:, 2:]
Y = data.iloc[:, 1].to_frame()

# 对特征进行标准化
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# PCA降维
pca = PCA()
plot_X = pca.fit_transform(X_scaled)

# 绘制所有数据点
fig, ax = plt.subplots()
norm = Normalize(vmin=Y.min(), vmax=Y.max())
scatter = ax.scatter(plot_X[:, 0], plot_X[:, 1], c=Y.values, cmap="coolwarm", norm=norm)
cb = fig.colorbar(scatter, extend='both')
cb.ax.set_title("Melting point(℃)")
# ax.set_title("")
ax.set_xlabel("PCA1")
ax.set_ylabel("PCA2")




# 绘制PCA1=-0.4的虚线
ax.axvline(x=-0.6, linestyle='--', color='k')
# 在屏幕上打印符合条件的样本点的数值
print("图一中 PCA1 < -0.4 的样本点的值:")
for i in range(len(plot_X)):
    if plot_X[i, 0] < -0.4:
        print(f"点 ({plot_X[i, 0]}, {plot_X[i, 1]}) 的值: {Y.values[i, 0]}")
        
        

# 将 PCA1 = -0.6 转换为原始特征空间中的值
pca1_value = -0.6
pca_components = pca.components_[0]  # 获取第一个主成分的方向
original_feature_values = pca1_value * pca_components * np.sqrt(pca.explained_variance_)

# 计算原始特征的方程
original_equation = " + ".join([f"{original_feature_values[i]:.3f} * Feature{i+1}" for i in range(len(original_feature_values))])

# 打印原始特征的方程
print("PCA1 = -0.6 对应的原始特征方程:")
print(original_equation)



# =============================================================================
# # 在点旁边标注值
# for i in range(len(plot_X)):
#     ax.text(plot_X[i, 0], plot_X[i, 1], str(Y.values[i, 0]), fontsize=8)
# =============================================================================

# 添加虚线框
ax.plot([-3.55, -2.6], [-0.39, -0.39], 'k--')
ax.plot([-3.55, -2.6], [0.59, 0.59], 'k--')
ax.plot([-3.55, -3.55], [-0.39, 0.59], 'k--')
ax.plot([-2.6, -2.6], [-0.39, 0.59], 'k--')

# 筛选出Y小于30的点
mask = Y.values.flatten() < 20
X_filtered = X[mask]
Y_filtered = Y[mask]

# 对筛选后的数据进行标准化
X_filtered_scaled = scaler.transform(X_filtered)

# 使用原始PCA对象对筛选后的数据进行降维
plot_X_filtered = pca.transform(X_filtered_scaled)

# 绘制筛选后的数据点
fig_filtered, ax_filtered = plt.subplots()
norm_filtered = Normalize(vmin=Y_filtered.min(), vmax=Y_filtered.max())
scatter_filtered = ax_filtered.scatter(plot_X_filtered[:, 0], plot_X_filtered[:, 1], c=Y_filtered.values, cmap="coolwarm", norm=norm_filtered)
cb_filtered = fig_filtered.colorbar(scatter_filtered, extend='both')
cb_filtered.ax.set_title("Melting point(℃)")
# ax_filtered.set_title("Y小于30的点的数据分布")
ax_filtered.set_xlabel("PCA1")
ax_filtered.set_ylabel("PCA2")

# 在筛选后的点旁边标注值
for i in range(len(plot_X_filtered)):
    ax_filtered.text(plot_X_filtered[i, 0], plot_X_filtered[i, 1], str(Y_filtered.values[i, 0]), fontsize=8)

# 添加虚线框
ax_filtered.plot([-3.55, -2.6], [-0.39, -0.39], 'k--')
ax_filtered.plot([-3.55, -2.6], [0.59, 0.59], 'k--')
ax_filtered.plot([-3.55, -3.55], [-0.39, 0.59], 'k--')
ax_filtered.plot([-2.6, -2.6], [-0.39, 0.59], 'k--')

# 打印框内的点的值
print("原始数据分布图中框内的点的值:")
for i in range(len(plot_X)):
    if -3.55 < plot_X[i, 0] < -2.6 and -0.39 < plot_X[i, 1] < 0.59:
        print(f"点 ({plot_X[i, 0]}, {plot_X[i, 1]}) 的值: {Y.values[i, 0]}")

print("\n筛选后数据分布图中框内的点的值:")
for i in range(len(plot_X_filtered)):
    if -3.55 < plot_X_filtered[i, 0] < -2.6 and -0.39 < plot_X_filtered[i, 1] < 0.59:
        print(f"点 ({plot_X_filtered[i, 0]}, {plot_X_filtered[i, 1]}) 的值: {Y_filtered.values[i, 0]}")

# 读取新数据集
new_data = pd.read_csv(r'C:/Users/senji/Desktop/模式识别ANN/实验验证样本集_填充描述符后ann模式识别实验样本.csv')
new_data.index = new_data.index.astype(str)

# 提取新数据集的特征
X_new = new_data.iloc[:, 2:]

# 对新数据集的特征进行标准化
X_new_scaled = scaler.transform(X_new)

# 使用原始PCA对象对新数据集进行降维
plot_X_new = pca.transform(X_new_scaled)

# 在原始数据分布图中绘制新数据点（五角星）
ax.scatter(plot_X_new[:, 0], plot_X_new[:, 1], color='black', marker='*', label='New Data')
ax.legend()

# 在筛选后的数据分布图中绘制新数据点（五角星）
ax_filtered.scatter(plot_X_new[:, 0], plot_X_new[:, 1], color='black', marker='*', label='New Data')
ax_filtered.legend()

plt.show()
#  -3.5<PCA1<-2.6,-0.39<PCA2<0.59