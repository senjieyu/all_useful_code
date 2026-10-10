# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 11:21:52 2023

@author: senji
"""

import pandas as pd
from sklearn.model_selection import train_test_split

# 读取原始数据文件
data = pd.read_csv("用于模型总数据.csv")

# 划分训练集和测试集
train_data, test_data = train_test_split(data, test_size=0.1, random_state=42)

# 保存划分后的数据为CSV文件
train_data.to_csv("train.csv", index=False)
test_data.to_csv("test.csv", index=False)

print("训练集和测试集已保存为 train.csv 和 test.csv")
