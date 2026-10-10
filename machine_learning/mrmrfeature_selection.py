# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 14:38:43 2023

@author: senji
"""

#数据格式，第一列为target  其余列为特征
#参数一：DataFrame，要求第一列是目标量（标注），其他列是特征量，其中首行必须是特征名称，且必须是字符形式，例如a1,a2，等
#参数二：选择的方法，有'MID'、'MIQ'两种。MID是基于互信息的mRMR，MIQ是基于商的mRMR。
#参数三：要求int类型，最后输出的特征数量。
import pandas as pd
import pymrmr

# Read in the CSV data file
data = pd.read_csv("8.18删去不可变特征_train.csv")

# Perform feature selection
mr = pymrmr.mRMR(data, 'MIQ', 20)

print(mr)

