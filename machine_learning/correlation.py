import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.sans-serif'] = ['SimHei']  # 设置字体为黑体

# 读取CSV文件
def read_csv(file_path):
    return pd.read_csv(file_path)

# 计算皮尔逊相关性
def calculate_pearson_correlation(df):
    correlations = {}
    features = df.columns[3:]  # 所有特征列
    target_column = df.columns[1]  # 目标列

    for feature in features:
        correlation = df[feature].corr(df[target_column], method='pearson')
        correlations[feature] = correlation

    return correlations

# 绘制柱状图
def plot_correlations(correlations):
    plt.figure(figsize=(10, 6), dpi=300)
    ax = plt.gca()
    ax.bar(correlations.keys(), correlations.values(), color='skyblue')
    plt.title('目标值与各特征之间的相关系数', fontproperties='SimHei',fontsize = 16)
    plt.xlabel('特征', fontproperties='SimHei',fontsize = 16)
    plt.ylabel('皮尔逊相关系数', fontproperties='SimHei',fontsize = 16)

    num_x_ticks = 9
    x_tick_indices = [int(i * (len(correlations) - 1) / (num_x_ticks - 1)) for i in range(num_x_ticks)]
    x_tick_labels = list(correlations.keys())[::len(correlations) // num_x_ticks]
    ax.set_xticks(x_tick_indices)
    ax.set_xticklabels(x_tick_labels, rotation=45, ha='right')  # 旋转45度，并水平对齐到右边

    plt.grid(axis='y')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # 输入CSV文件路径
    file_path = r"C:/Users/senji/Desktop/梅菜扣肉/一阶导平滑.csv"

    # 读取CSV文件
    df = read_csv(file_path)

    # 计算皮尔逊相关性
    correlations = calculate_pearson_correlation(df)

    # 绘制柱状图
    plot_correlations(correlations)
