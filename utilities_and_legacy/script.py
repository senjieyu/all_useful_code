# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 17:11:13 2023

@author: senji
"""

import pandas as pd

# 读取原始CSV文件，假设文件名为 'data.csv'
data = pd.read_csv('test.csv')

# 指定要提取的特征列名（假设要提取 'feature1' 和 'feature3'）
selected_feature_names = ['No','cop_dataprocess','RealtimeLoad','CWTempReturn01','CWTempReturnEnterChiller01','SmallChWReturnTemp','ChWTempSupply01','CWTempSupply01','OutdoorDewPoinTemp']  # 修改为你需要提取的特征列名

# 根据选定的特征列名提取特征数据
selected_features = data[selected_feature_names]

# 保存提取的特征数据到新的CSV文件
selected_features.to_csv('8.18mrmr前20，RF前7_test.csv', index=False)
