# -*- coding: utf-8 -*-
"""
优化后的SHAP可视化代码（字体统一调整版）
"""
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
import pandas as pd
import shap
import matplotlib.pyplot as plt
import numpy as np

# ==================== 全局样式设置 ====================
plt.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 24,                    # 基础字体大小统一设为24
    "axes.titlesize": 24,               # 主标题大小
    "axes.labelsize": 24,               # 坐标轴标签大小
    "xtick.labelsize": 24,              # X轴刻度
    "ytick.labelsize": 24,              # Y轴刻度
    "figure.dpi": 1000,                 # 超高分辨率
    "xtick.direction": "in",            # 刻度朝内
    "ytick.direction": "in",
    "axes.titlepad": 20                 # 标题间距
})

# ==================== 数据准备 ====================
train = pd.read_csv("train副本.csv").iloc[:, 1:].astype(float)
test = pd.read_csv("test副本.csv").iloc[:, 1:].astype(float)

X_train = train.drop(columns=["Melting_point"])
y_train = train["Melting_point"]
X_test = test.drop(columns=["Melting_point"])
y_test = test["Melting_point"]

# ==================== 模型训练 ====================
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

model = SVR(C=1000.0, gamma=0.12, epsilon=0.01, kernel="rbf")
model.fit(X_train_scaled, y_train)

# ==================== SHAP可视化 ====================
background = shap.kmeans(X_train_scaled, 10)
explainer = shap.KernelExplainer(model.predict, background)
shap_values_train = explainer.shap_values(X_train_scaled, silent=True)

# 1. SHAP总结图（点图）
# 1. SHAP总结图（点图）
plt.figure(figsize=(12, 8))
shap.summary_plot(shap_values_train, features=X_train, 
                  feature_names=X_train.columns,
                  show=False)

# 设置主图字体
ax = plt.gca()
for item in ax.findobj(match=plt.Text):
    item.set_fontsize(24)
    item.set_fontname("Times New Roman")

# 设置 colorbar 字体
fig = plt.gcf()
for ax in fig.axes:
    if hasattr(ax, 'get_ylabel') and 'Feature value' in ax.get_ylabel():
        ax.tick_params(labelsize=24)
        ax.yaxis.label.set_size(24)
        ax.yaxis.label.set_fontname("Times New Roman")

plt.tight_layout()
plt.savefig("shap_summary_dot.tif", dpi=1000, bbox_inches='tight')
plt.show()


# 2. SHAP条形图（重点调整部分）
plt.figure(figsize=(12, 8))
shap.summary_plot(shap_values_train, features=X_train, 
                 feature_names=X_train.columns,
                 plot_type="bar",
                 show=False)

ax = plt.gca()

# 统一设置所有文本（包括条形标签和坐标轴）
for item in ax.findobj(match=plt.Text):
    item.set_fontsize(24)
    item.set_fontname("Times New Roman")

# 特别确保X轴标题
ax.set_xlabel("Average SHAP Value Magnitude", fontsize=24)

# 设置条形颜色和边框
for patch in ax.patches:
    patch.set_facecolor((245/255, 180/255, 130/255))
    patch.set_edgecolor('black')
    patch.set_linewidth(0.5)

plt.tight_layout()
plt.show()