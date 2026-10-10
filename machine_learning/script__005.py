import pandas as pd
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, median_absolute_error
from scipy.stats import pearsonr
from autogluon.tabular import TabularPredictor

# 设置 Pandas 显示选项
pd.set_option('display.max_columns', None)  # 显示所有列
pd.set_option('display.max_rows', None)     # 显示所有行
pd.set_option('display.width', None)        # 自动调整宽度
pd.set_option('display.max_colwidth', None) # 显示每列的完整内容

# 读取数据
data = pd.read_csv('C:/Users/senji/Desktop/测试自动机器学习.csv')

# 提取特征和标签
label_column = data.columns[1]  # 标签列

# 划分数据集
train_data = data.sample(frac=0.8, random_state=42)
test_data = data.drop(train_data.index)

# 训练回归模型
predictor = TabularPredictor(label=label_column).fit(train_data)

# 输出领导者板
leaderboard = predictor.leaderboard(silent=True)
print(leaderboard)

# 进行预测
train_predictions = predictor.predict(train_data)
test_predictions = predictor.predict(test_data)

# 计算训练集性能评估指标
train_rmse = mean_squared_error(train_data[label_column], train_predictions, squared=False)  # RMSE
train_mse = mean_squared_error(train_data[label_column], train_predictions)  # MSE
train_mae = mean_absolute_error(train_data[label_column], train_predictions)  # MAE
train_r2 = r2_score(train_data[label_column], train_predictions)  # R²
train_median_ae = median_absolute_error(train_data[label_column], train_predictions)  # 中位数绝对误差
train_pearson_corr, _ = pearsonr(train_data[label_column], train_predictions)  # 皮尔逊相关系数

# 计算测试集性能评估指标
test_rmse = mean_squared_error(test_data[label_column], test_predictions, squared=False)  # RMSE
test_mse = mean_squared_error(test_data[label_column], test_predictions)  # MSE
test_mae = mean_absolute_error(test_data[label_column], test_predictions)  # MAE
test_r2 = r2_score(test_data[label_column], test_predictions)  # R²
test_median_ae = median_absolute_error(test_data[label_column], test_predictions)  # 中位数绝对误差
test_pearson_corr, _ = pearsonr(test_data[label_column], test_predictions)  # 皮尔逊相关系数

# 打印训练集和测试集结果
print("训练集性能评估指标:")
print(f"均方根误差 (RMSE): {train_rmse}")
print(f"均方误差 (MSE): {train_mse}")
print(f"平均绝对误差 (MAE): {train_mae}")
print(f"决定系数 (R²): {train_r2}")
print(f"中位数绝对误差: {train_median_ae}")
print(f"皮尔逊相关系数: {train_pearson_corr}")

print("\n测试集性能评估指标:")
print(f"均方根误差 (RMSE): {test_rmse}")
print(f"均方误差 (MSE): {test_mse}")
print(f"平均绝对误差 (MAE): {test_mae}")
print(f"决定系数 (R²): {test_r2}")
print(f"中位数绝对误差: {test_median_ae}")
print(f"皮尔逊相关系数: {test_pearson_corr}")
