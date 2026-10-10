# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 15:12:19 2024

@author: senji
"""

# =============================================================================
# # import pandas as pd
# # import pymrmr
# # # 第一列是每个样本的分类（目标）变量。其余列是算法可以选择的不同变量（特征）。
# 
# # # Read in the CSV data file
# # data = pd.read_excel(r"D:/公司/卢天-纪晓波25-3-17/纪晓波20250216/分类问题/特征筛选mrmr.xlsx")
# 
# # # Perform feature selection
# # mr = pymrmr.mRMR(data, 'MIQ', 40)
# 
# # print(mr)
# =============================================================================

import pandas as pd
import pymrmr

# 读取 Excel 文件
file_path = r"D:/公司/卢天-纪晓波25-3-17/纪晓波20250216/分类问题/特征筛选mrmr.xlsx"
data = pd.read_excel(file_path)

# 删除常数列
constant_columns = data.columns[data.nunique() == 1]  # 找到常数列
data = data.drop(columns=constant_columns)  # 删除常数列

# 检查删除常数列后的数据
print("删除常数列后的数据：")
print(data.head())

# 检查数据类型
print("数据类型检查：")
print(data.dtypes)

# 确保所有特征列都是数值类型
non_numeric_columns = data.select_dtypes(exclude=['number']).columns
if len(non_numeric_columns) > 0:
    print("发现非数值列，正在删除：", non_numeric_columns)
    data = data.drop(columns=non_numeric_columns)

# 确保第三个参数是整数
num_features = int(40)  # 显式转换为整数

# 检查数据类型
print("数据类型检查：")
print("data 类型:", type(data))  # 应该是 <class 'pandas.core.frame.DataFrame'>
print("num_features 类型:", type(num_features))  # 应该是 <class 'int'>

# 执行 mRMR 特征选择
try:
    selected_features = pymrmr.mRMR(data, 'MIQ', num_features)
    print("选中的特征：", selected_features)
except Exception as e:
    print("执行 mRMR 时出错:", e)

# 将选中的特征列提取出来
if 'selected_features' in locals():
    selected_data = data[selected_features]

    # 将结果保存到 Excel 文件
    output_path = r"D:/公司/卢天-纪晓波25-3-17/纪晓波20250216/分类问题/特征筛选mrmr_selected.xlsx"
    selected_data.to_excel(output_path, index=False)
    print(f"选中的特征已保存到: {output_path}")