import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set Times New Roman as the default font
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.sans-serif'] = ['SimHei']  # 设置字体为黑体

# 读取CSV文件，假设第一列为索引
data = pd.read_csv('C:/Users/senji/Desktop/梅菜扣肉/一阶导函数.csv', index_col=0)

# 获取数据的行数和列数
num_rows, num_cols = data.shape

# 创建一个图形和轴对象，设置dpi为300
fig, ax = plt.subplots(dpi=300)

# 对每一行数据进行作图
for i in range(num_rows):
    row_data = data.iloc[i]  # 获取第i行数据
    ax.plot(row_data, label=f'Sample {i+1}')  # 绘制折线图，label为样本的编号

# 设置x轴刻度
num_x_ticks = 8
x_tick_indices = [int(i * (num_cols - 1) / (num_x_ticks - 1)) for i in range(num_x_ticks)]
x_tick_labels = [data.columns[idx] for idx in x_tick_indices]
ax.set_xticks(x_tick_indices)
ax.set_xticklabels(x_tick_labels, rotation=45)  # 旋转45度

# 添加标题和标签，设置字体为SimHei
ax.set_title('Line Plot for Each Sample', fontproperties='SimHei')
ax.set_xlabel('波长', fontproperties='SimHei')  # 设置x轴标签
ax.set_ylabel('强度', fontproperties='SimHei')  # 设置y轴标签

# 显示图例
# ax.legend()

# 显示图形
plt.show()
