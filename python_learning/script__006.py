import os
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from tqdm import tqdm  # 引入tqdm库，用于显示进度条


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


# 定义卷积神经网络模型
class EEGNet(nn.Module):
    def __init__(self):
        super(EEGNet, self).__init__()
        self.conv1 = nn.Conv2d(1, 32, kernel_size=(3, 3), padding=1)
        self.pool = nn.MaxPool2d(kernel_size=(2, 2), stride=2, padding=0)
        self.conv2 = nn.Conv2d(32, 64, kernel_size=(3, 3), padding=1)

        # 根据池化操作，动态计算全连接层的输入大小
        self._to_linear = None
        self.convs(torch.randn(1, 1, 2500, 32))  # 通过一个假输入来计算self._to_linear

        self.fc1 = nn.Linear(self._to_linear, 128)
        self.fc2 = nn.Linear(128, 2)

    def convs(self, x):
        x = self.pool(torch.relu(self.conv1(x)))
        x = self.pool(torch.relu(self.conv2(x)))
        if self._to_linear is None:
            self._to_linear = int(np.prod(x.size()[1:]))
        return x

    def forward(self, x):
        x = self.convs(x)
        x = x.view(x.size(0), -1)  # 展平
        x = torch.relu(self.fc1(x))
        x = self.fc2(x)
        return x


# 设置路径
base_path = r'D:\EEG\后2500行，分类数据集'

# 载入数据
data, labels = load_data(base_path)
# 添加一个通道维度
data = data[:, None, :, :]  # shape (num_samples, 1, num_rows, num_cols)
X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=0.2, random_state=42)

train_dataset = EEGDataset(X_train, y_train)
test_dataset = EEGDataset(X_test, y_test)

train_loader = DataLoader(train_dataset, batch_size=16, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=16, shuffle=False)

# 初始化模型、损失函数和优化器
model = EEGNet()
criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# 训练模型
num_epochs = 20
for epoch in range(num_epochs):
    model.train()
    running_loss = 0.0
    correct = 0
    total = 0
    for inputs, labels in tqdm(train_loader, desc=f"Epoch {epoch + 1}/{num_epochs}", unit="batch"):
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
        running_loss += loss.item()

        # Calculate training accuracy
        _, preds = torch.max(outputs, 1)
        correct += (preds == labels).sum().item()
        total += labels.size(0)

    train_accuracy = correct / total
    print(f'Epoch [{epoch + 1}/{num_epochs}], Loss: {running_loss / len(train_loader):.4f}, Train Accuracy: {train_accuracy:.4f}')

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
