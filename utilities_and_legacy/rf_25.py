# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 14:50:32 2023

@author: senji
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 17:30:42 2023

@author: senji
"""

import pandas as pd
from sklearn.ensemble import RandomForestRegressor

# 读取CSV文件，假设文件名为 'data.csv'
data = pd.read_csv('8.18前20_train.csv')

# 提取目标变量（第二列）
target = data.iloc[:, 1]

# 提取特征变量（第三列到最后一列）
features = data.iloc[:, 2:]

# 创建随机森林回归器
rf_regressor = RandomForestRegressor(n_estimators=100, random_state=42)

# 训练随机森林回归模型
rf_regressor.fit(features, target)

# 获取特征重要性评分
feature_importances = rf_regressor.feature_importances_

# 将特征重要性评分和对应的特征列名组合成一个DataFrame
feature_importance_df = pd.DataFrame({'Feature': features.columns, 'Importance': feature_importances})

# 按照重要性降序排列
feature_importance_df = feature_importance_df.sort_values(by='Importance', ascending=False)

# 输出所有特征的重要性得分
print("Feature Importance Scores:")
print(feature_importance_df)

# 将选定的特征和目标变量重新整合成一个新的DataFrame
#selected_data = data[['No']+['cop_dataprocess'] + list(feature_importance_df['Feature'])]

# 保存选定的特征和目标变量到新的CSV文件
#selected_data.to_csv('.csv', index=False)
