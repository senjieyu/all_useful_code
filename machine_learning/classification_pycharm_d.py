#  '运行在虚拟环境上，虚拟环境配置地址：D:\EEG\EEG_pycharm_code\venv\Scripts\Python.exe'
import pandas as pd
from lazypredict.Supervised import LazyClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

# 读取数据集
data = pd.read_csv(r'D:/脑电/平均值数据集.csv')

# 提取特征和目标变量
X, y = data.iloc[:, 1:], data.iloc[:, 0]

# 分割数据集为训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 初始化 LazyClassifier
clf = LazyClassifier(predictions=True)

# 拟合模型
models, predictions = clf.fit(X_train, X_test, y_train, y_test)

# 打印每种模型的训练集和测试集分类准确率
print("Model Performance:")
print(models.to_string())


# 将结果保存到一个csv文件
results_df = models.copy()
results_df['Test Accuracy'] = [accuracy_score(y_test, predictions[model_name]) for model_name in models.index]
results_df.to_csv('model_performance.csv', index=True)
