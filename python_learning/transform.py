import os
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from tqdm import tqdm

# 定义数据集类
class EEGDataset(Dataset):
    def __init__(self, data, labels):
        self.data = data
        self.labels = labels

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        sample = self.data[idx]
        return torch.tensor(sample, dtype=torch.float32), torch.tensor(self.labels[idx], dtype=torch.long)

# 读取数据并创建数据集
def load_data(base_path):
    data = []
    labels = []
    for label in ['1', '2']:
        folder = os.path.join(base_path, label)
        for file_name in os.listdir(folder):
            if file_name.endswith('.txt'):
                file_path = os.path.join(folder, file_name)
                file_data = np.loadtxt(file_path, delimiter=',', usecols=range(1, 33))  # 读取第2-33列
                data.append(file_data)  # 将文件中的所有数据作为一个二维数组
                labels.append(int(label) - 1)  # 标签从0开始
    return np.array(data), np.array(labels)

# 定义Transformer模型
class EEGTransformer(nn.Module):
    def __init__(self, input_dim, num_heads, num_layers, num_classes):
        super(EEGTransformer, self).__init__()
        self.embedding = nn.Linear(input_dim, 512)
        encoder_layer = nn.TransformerEncoderLayer(d_model=512, nhead=num_heads, batch_first=True)
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        self.fc = nn.Linear(512, num_classes)

    def forward(self, x):
        x = self.embedding(x)
        x = self.transformer(x)
        x = x.mean(dim=1)  # 取所有时间步的平均值
        x = self.fc(x)
        return x

# 设置路径
base_path = r'D:\EEG\后2500行，分类数据集'

# 载入数据
data, labels = load_data(base_path)
X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=0.2, random_state=42)

train_dataset = EEGDataset(X_train, y_train)
test_dataset = EEGDataset(X_test, y_test)

train_loader = DataLoader(train_dataset, batch_size=16, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=16, shuffle=False)

# 初始化模型、损失函数和优化器
input_dim = 32
num_heads = 8
num_layers = 2
num_classes = 2
model = EEGTransformer(input_dim, num_heads, num_layers, num_classes)

criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# 训练模型
num_epochs = 20
for epoch in range(num_epochs):
    model.train()
    running_loss = 0.0
    for inputs, labels in tqdm(train_loader, desc=f'Epoch {epoch + 1}/{num_epochs}', unit='batch'):
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
        running_loss += loss.item()

    print(f'Epoch [{epoch + 1}/{num_epochs}], Loss: {running_loss / len(train_loader):.4f}')

    # 打印训练集的准确率
    model.eval()
    train_preds = []
    train_labels = []
    with torch.no_grad():
        for inputs, labels in train_loader:
            outputs = model(inputs)
            _, preds = torch.max(outputs, 1)
            train_preds.extend(preds.numpy())
            train_labels.extend(labels.numpy())
    train_accuracy = accuracy_score(train_labels, train_preds)
    print(f'Training Accuracy: {train_accuracy:.4f}')

# 测试模型
model.eval()
all_preds = []
all_labels = []
with torch.no_grad():
    for inputs, labels in test_loader:
        outputs = model(inputs)
        _, preds = torch.max(outputs, 1)
        all_preds.extend(preds.numpy())
        all_labels.extend(labels.numpy())

accuracy = accuracy_score(all_labels, all_preds)
print(f'Test Accuracy: {accuracy:.4f}')
