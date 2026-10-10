import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.cross_decomposition import PLSRegression
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import accuracy_score

from matplotlib import rcParams

# Set Times New Roman as the default font
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']

# 读取CSV文件



df = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/部位1_不含重现样本.csv')


# 提取目标值列和特征列
X = df.iloc[:, 2:]  # 假设特征列从第三列开始
y = df.iloc[:, 1]   # 假设目标值列为第二列

# 初始化最佳结果变量
results = {
    "SVM": {"dimensions": [], "train_accuracy": [], "test_accuracy": []},
    "Logistic Regression": {"dimensions": [], "train_accuracy": [], "test_accuracy": []},
    "K-Nearest Neighbors": {"dimensions": [], "train_accuracy": [], "test_accuracy": []},
    "Decision Tree": {"dimensions": [], "train_accuracy": [], "test_accuracy": []},
    "Random Forest": {"dimensions": [], "train_accuracy": [], "test_accuracy": []},
    "Gradient Boosting": {"dimensions": [], "train_accuracy": [], "test_accuracy": []}
}

classifiers = {
    "SVM": SVC(),
    "Logistic Regression": LogisticRegression(max_iter=1000),
    "K-Nearest Neighbors": KNeighborsClassifier(),
    "Decision Tree": DecisionTreeClassifier(random_state=42),
    "Random Forest": RandomForestClassifier(random_state=42),
    "Gradient Boosting": GradientBoostingClassifier(random_state=42)
}

for classifier_name, classifier in classifiers.items():
    for n_components in range(2, 25):  # 尝试2到100个维度
        pls = PLSRegression(n_components=n_components)

        # 划分训练集和测试集
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        # 训练PLS模型
        pls.fit(X_train, y_train)

        # 使用训练好的模型对数据进行降维
        X_train_pls = pls.transform(X_train)
        X_test_pls = pls.transform(X_test)

        # 训练分类器模型
        classifier.fit(X_train_pls, y_train)

        # 预测测试集
        y_pred = classifier.predict(X_test_pls)

        # 计算测试集准确率
        test_accuracy = accuracy_score(y_test, y_pred)

        # 计算训练集准确率
        train_accuracy = accuracy_score(y_train, classifier.predict(X_train_pls))

        # 保存结果
        results[classifier_name]["dimensions"].append(n_components)
        results[classifier_name]["train_accuracy"].append(train_accuracy)
        results[classifier_name]["test_accuracy"].append(test_accuracy)

# =============================================================================
# # 绘制点线图
# plt.figure(figsize=(10, 6),dpi = 300)
# 
# 
# colors = ['#8ECFC9', '#FFBE7A', '#FA7F6F', '#82B0D2', '#BEB8DC', '#E7DAD2']
# color_index = 0
# 
# for classifier_name, result in results.items():
#     plt.plot(result["dimensions"], result["train_accuracy"], marker='o', color=colors[color_index], label=f'{classifier_name} Train Accuracy')
#     color_index = (color_index + 1) % len(colors)  # 循环使用颜色列表中的颜色
# plt.legend()
# plt.xlabel('Dimensions',fontsize = 24,)
# plt.ylabel('Accuracy',fontsize = 24,)
# plt.title('Train Accuracy',fontsize = 24,)
# plt.tick_params(axis='both', labelsize=16, )
# plt.show()
# 
# plt.figure(figsize=(10, 6),dpi = 300)
# for classifier_name, result in results.items():
#     plt.plot(result["dimensions"], result["test_accuracy"], marker='o', linestyle='--', color=colors[color_index], label=f'{classifier_name} Test Accuracy')
#     color_index = (color_index + 1) % len(colors)  # 循环使用颜色列表中的颜色
# 
# plt.legend()
# plt.xlabel('Dimensions',fontsize = 24,)
# plt.ylabel('Accuracy',fontsize = 24,)
# plt.title('Test Accuracy',fontsize = 24,)
# plt.tick_params(axis='both', labelsize=16, )
# plt.show()
# 
# =============================================================================

import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(18, 6), dpi=300)
colors = ['#8ECFC9', '#FFBE7A', '#FA7F6F', '#82B0D2', '#BEB8DC', '#E7DAD2']
# 绘制训练准确率曲线
axes[0].set_prop_cycle('color', colors)
for classifier_name, result in results.items():
    axes[0].plot(result["dimensions"], result["train_accuracy"], marker='o', label=f'{classifier_name} Train Accuracy')
axes[0].legend(fontsize=14)
axes[0].set_xlabel('Dimensions', fontsize=24)
axes[0].set_ylabel('Accuracy', fontsize=24)
axes[0].set_title('Train Accuracy', fontsize=24)
axes[0].tick_params(axis='both', labelsize=16)

# 绘制测试准确率曲线
axes[1].set_prop_cycle('color', colors)
for classifier_name, result in results.items():
    axes[1].plot(result["dimensions"], result["test_accuracy"], marker='o', linestyle='--', label=f'{classifier_name} Test Accuracy')
axes[1].legend(fontsize=14)
axes[1].set_xlabel('Dimensions', fontsize=24)
axes[1].set_ylabel('Accuracy', fontsize=24)
axes[1].set_title('Test Accuracy', fontsize=24)
axes[1].tick_params(axis='both', labelsize=16)

plt.tight_layout()  # 自动调整子图布局，避免重叠
plt.show()







