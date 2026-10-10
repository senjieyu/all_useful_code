import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
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
df = pd.read_csv(r'C:/Users/senji/Desktop/烟草仅部位/一阶导函数.csv')

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
    "SVM": SVC(random_state=42),
    "Logistic Regression": LogisticRegression(max_iter=1000, random_state=42),
    "K-Nearest Neighbors": KNeighborsClassifier(),
    "Decision Tree": DecisionTreeClassifier(random_state=42),
    "Random Forest": RandomForestClassifier(random_state=42),
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
best_X_train_pca = None
best_X_test_pca = None

for classifier_name, classifier in classifiers.items():
    for n_components in range(2, 25):  # 尝试2到100个维度
        pca = PCA(n_components=n_components)  # 使用PCA进行降维

        # 划分训练集和测试集
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        # 使用PCA对数据进行降维
        X_train_pca = pca.fit_transform(X_train)
        X_test_pca = pca.transform(X_test)

        # 训练分类器模型
        classifier.fit(X_train_pca, y_train)

        # 预测测试集
        y_pred = classifier.predict(X_test_pca)

        # 计算测试集准确率
        test_accuracy = accuracy_score(y_test, y_pred)

        # 计算训练集准确率
        train_accuracy = accuracy_score(y_train, classifier.predict(X_train_pca))

        # 保存结果
        results[classifier_name]["dimensions"].append(n_components)
        results[classifier_name]["train_accuracy"].append(train_accuracy * 100)  # 将准确率转换为百分比形式
        results[classifier_name]["test_accuracy"].append(test_accuracy * 100)  # 将准确率转换为百分比形式

        # 输出结果
        print(f"Classifier: {classifier_name}, Dimensions: {n_components}, Train Accuracy: {train_accuracy * 100:.2f}%, Test Accuracy: {test_accuracy * 100:.2f}%")

        # 保存最佳结果
        if test_accuracy > max_accuracy:
            max_accuracy = test_accuracy
            best_classifier = classifier_name
            best_dimensions = n_components
            best_X_train = X_train_pca
            best_X_test = X_test_pca
            best_y_train = y_train
            best_y_test = y_test
            best_y_pred = y_pred
            best_X_test_pca = X_test_pca

print(f"Best Dimensionality Reduction Result - Dimensions: {best_dimensions}, Classifier: {best_classifier}, 准确率: {max_accuracy * 100:.2f}%")

# 绘制折线图
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

# 绘制混淆矩阵
conf_matrix = confusion_matrix(best_y_test, best_y_pred)
plt.figure(figsize=(8, 6))
plt.imshow(conf_matrix, interpolation='nearest', cmap=plt.cm.Blues)
plt.title('Confusion Matrix')
plt.colorbar()
tick_marks = range(len(set(best_y_test)))
plt.xticks(tick_marks, set(best_y_test), rotation=45)
plt.yticks(tick_marks, set(best_y_test))
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.show()

# 打印混淆矩阵
print("测试集混淆矩阵:")
print(conf_matrix)
