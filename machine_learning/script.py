import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.colors import Normalize

plt.rcParams["figure.dpi"] = 300
plt.rcParams["font.sans-serif"] = ["SimHei"]
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
cb.ax.set_title("COP")
ax.set_title("删除融冰/蓄冰后的数据分布")
ax.set_xlabel("PCA1")
ax.set_ylabel("PCA2")

# 在点旁边标注值
for i in range(len(plot_X)):
    ax.text(plot_X[i, 0], plot_X[i, 1], str(Y.values[i, 0]), fontsize=8)

# 添加虚线框
ax.plot([-3.5, -2.6], [-0.39, -0.39], 'k--')
ax.plot([-3.5, -2.6], [0.59, 0.59], 'k--')
ax.plot([-3.5, -3.5], [-0.39, 0.59], 'k--')
ax.plot([-2.6, -2.6], [-0.39, 0.59], 'k--')

# 筛选出Y小于20的点
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
cb_filtered.ax.set_title("COP")
ax_filtered.set_title("Y小于20的点的数据分布")
ax_filtered.set_xlabel("PCA1")
ax_filtered.set_ylabel("PCA2")

# 在筛选后的点旁边标注值
for i in range(len(plot_X_filtered)):
    ax_filtered.text(plot_X_filtered[i, 0], plot_X_filtered[i, 1], str(Y_filtered.values[i, 0]), fontsize=8)

# 添加虚线框
ax_filtered.plot([-3.5, -2.6], [-0.39, -0.39], 'k--')
ax_filtered.plot([-3.5, -2.6], [0.59, 0.59], 'k--')
ax_filtered.plot([-3.5, -3.5], [-0.39, 0.59], 'k--')
ax_filtered.plot([-2.6, -2.6], [-0.39, 0.59], 'k--')

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
