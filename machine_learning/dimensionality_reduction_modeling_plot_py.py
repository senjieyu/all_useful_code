import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.cross_decomposition import PLSRegression
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import accuracy_score, confusion_matrix
from matplotlib import rcParams

# Set Times New Roman as the default font
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']

# 读取CSV文件
df = pd.read_csv(r'D:/EEG/标准差数据集.csv')

# 提取目标值列和特征列
X = df.iloc[:, 1:]  # 假设特征列从第三列开始
y = df.iloc[:, 0]   # 假设目标值列为第二列

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
    "SVM": SVC(random_state=42),
    "Logistic Regression": LogisticRegression(max_iter=2000,random_state=42),
    "K-Nearest Neighbors": KNeighborsClassifier(),
    "Decision Tree": DecisionTreeClassifier(random_state=42),
    "Random Forest": RandomForestClassifier(random_state=3),#
    "Gradient Boosting": GradientBoostingClassifier(random_state=42)
}

max_accuracy = 0
best_classifier = None
best_dimensions = 0
best_X_train = None
best_X_test = None
best_y_train = None
best_y_test = None
best_y_pred = None
best_X_train_pls = None
best_X_test_pls = None

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
        results[classifier_name]["train_accuracy"].append(train_accuracy * 100)  # 将准确率转换为百分比形式
        results[classifier_name]["test_accuracy"].append(test_accuracy * 100)  # 将准确率转换为百分比形式
        if test_accuracy > max_accuracy:
            max_accuracy = test_accuracy
            best_classifier = classifier_name
            best_X_train_pls = X_train_pls
            best_dimensions = n_components
            best_X_train = X_train_pls
            best_X_test = X_test_pls
            best_y_train = y_train
            best_y_test = y_test
            best_y_pred = y_pred
            best_X_test_pls = X_test_pls


print(f"Best Dimensionality Reduction Result - Dimensions: {best_dimensions}, Classifier: {best_classifier}, 准确率：{max_accuracy * 100:.2f}%")

#111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111


import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(18, 6), dpi=300)
colors = ['#8ECFC9', '#FFBE7A', '#FA7F6F', '#82B0D2', '#BEB8DC', '#E7DAD2']
# 绘制训练准确率曲线
axes[0].set_prop_cycle('color', colors)
for classifier_name, result in results.items():
    axes[0].plot(result["dimensions"], result["train_accuracy"], marker='o', label=f'{classifier_name} Train Accuracy')
axes[0].legend(fontsize=14)
axes[0].set_xlabel('Dimensions', fontsize=24)
axes[0].set_ylabel('Accuracy (%)', fontsize=24)  # 修改 y 轴标签
axes[0].set_title('Train Accuracy', fontsize=24)
axes[0].tick_params(axis='both', labelsize=16)

# 绘制测试准确率曲线
axes[1].set_prop_cycle('color', colors)
for classifier_name, result in results.items():
    axes[1].plot(result["dimensions"], result["test_accuracy"], marker='o', linestyle='--', label=f'{classifier_name} Test Accuracy')
axes[1].legend(fontsize=14)
axes[1].set_xlabel('Dimensions', fontsize=24)
axes[1].set_ylabel('Accuracy (%)', fontsize=24)  # 修改 y 轴标签
axes[1].set_title('Test Accuracy', fontsize=24)
axes[1].tick_params(axis='both', labelsize=16)

plt.tight_layout()  # 自动调整子图布局，避免重叠
plt.show()

#111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111


best_X_train_pls_df = pd.DataFrame(best_X_train_pls, columns=[f"PLS_Component_{i+1}" for i in range(best_dimensions)])
best_X_test_pls_df = pd.DataFrame(best_X_test_pls, columns=[f"PLS_Component_{i+1}" for i in range(best_dimensions)])

# 合并特征数据和目标数据
combined_train_pls = pd.concat([best_y_train.reset_index(drop=True), best_X_train_pls_df], axis=1)
combined_test_pls = pd.concat([best_y_test.reset_index(drop=True), best_X_test_pls_df], axis=1)


#111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111


# 提取特征列和目标列
X_train = combined_train_pls.iloc[:, 1:]  # 特征列，从第二列开始到最后一列
y_train = combined_train_pls.iloc[:, 0]   # 目标列，第一列

X_test = combined_test_pls.iloc[:, 1:]
y_test = combined_test_pls.iloc[:, 0]


classifier = classifiers.get(best_classifier)


# 拟合模型
classifier.fit(X_train, y_train)

# 预测训练集和测试集
train_predictions = classifier.predict(X_train)
test_predictions = classifier.predict(X_test)

# 计算混淆矩阵
train_confusion_matrix = confusion_matrix(y_train, train_predictions)
test_confusion_matrix = confusion_matrix(y_test, test_predictions)

print("训练集混淆矩阵:")
print(train_confusion_matrix)
print("\n测试集混淆矩阵:")
print(test_confusion_matrix)


train_accuracy = accuracy_score(y_train, train_predictions) * 100  # 将训练集准确率转换为百分比形式
test_accuracy = accuracy_score(y_test, test_predictions) * 100  # 将测试集准确率转换为百分比形式

print("训练集准确率:", train_accuracy)
print("测试集准确率:", test_accuracy)


conf_matrix = confusion_matrix(best_y_test, best_y_pred)
plt.figure(figsize=(8, 6), dpi=300)
plt.imshow(conf_matrix, interpolation='nearest', cmap=plt.cm.Blues)
plt.title('Confusion Matrix', fontsize=24)
colorbar = plt.colorbar()
colorbar.ax.tick_params(labelsize=20)  # 调整右侧标尺的字体大小
tick_marks = range(len(set(best_y_test)))
plt.xticks(tick_marks, set(best_y_test), fontsize=20)
plt.yticks(tick_marks, set(best_y_test), fontsize=20)
plt.xlabel('Predicted Label', fontsize=24)
plt.ylabel('True Label', fontsize=24)

# 将混淆矩阵中的数字添加到每个方格的中心
for i in range(conf_matrix.shape[0]):
    for j in range(conf_matrix.shape[1]):
        plt.text(j, i, str(conf_matrix[i, j]), horizontalalignment="center", verticalalignment="center", fontsize=24)

plt.show()