import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
# 第一列数据为目标值列，其余列为特征列


# 读取CSV文件
df = pd.read_excel('C:/Users/senji/Desktop/工作簿1.xlsx').iloc[:50,:]  # 替换为你的文件路径


# 分离特征和目标变量
y = df.iloc[:, 0]  # 假设第一列是目标变量Bandgap
X = df.iloc[:, 1:]  # 其余列作为特征

# 标准化特征
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 应用PCA，计算第一主成分
pca = PCA(n_components=1)
X_pca = pca.fit_transform(X_scaled)

# 高斯过程回归
kernel = C(1.0, (1e-2, 1e2)) * RBF(10.0, (1e-1, 1e1))
gpr = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10, alpha=0.1)
gpr.fit(X_pca, y)

# 预测均值和标准差
x = np.linspace(X_pca.min(), X_pca.max(), 100).reshape(-1, 1)
y_pred, sigma = gpr.predict(x, return_std=True)

# 绘制散点图和拟合曲线
plt.figure(figsize=(10, 6), dpi=300)
scatter = plt.scatter(X_pca[:, 0], y, c=y, cmap='coolwarm', vmin=y.min(), vmax=y.max(), edgecolor='k', s=50)
plt.colorbar(scatter, label='Bandgap (eV)')

# 绘制高斯过程的预测均值和置信区间
plt.plot(x, y_pred, 'orange', lw=2, label='GPR Mean')
plt.fill_between(x[:, 0], y_pred - 1.96 * sigma, y_pred + 1.96 * sigma, color='orange', alpha=0.3, label='95% Confidence Interval')

plt.title('Proactive Searching Progress')
plt.xlabel('PCA 1')
plt.ylabel('Melting Point (℃)')
plt.legend()
# plt.grid(True)
plt.show()