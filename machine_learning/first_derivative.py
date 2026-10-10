# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 10:09:47 2024

@author: senji
"""

import pandas as pd
import numpy as np

# 读取原始数据集的CSV文件
df = pd.read_csv(r'C:/Users/senji/Desktop/梅菜扣肉/删除15异常值.csv')

# 提取光谱数据列（假设光谱数据从第二列开始）
spectra_data = df.iloc[:, 3:]

# 计算一阶导数
derivatives = np.gradient(spectra_data, axis=1)

# 构建包含导数数据的新DataFrame
derivatives_df = pd.DataFrame(derivatives, columns=spectra_data.columns)

# 将原始数据和导数数据合并为一个新的DataFrame
merged_df = pd.concat([df.iloc[:, :3], derivatives_df], axis=1)

# 保存导数数据到CSV文件
merged_df.to_csv('一阶导函数.csv', index=False)