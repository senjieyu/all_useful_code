import torch
from torch.utils.data import TensorDataset
a = torch.randn(500,2,3)
b = torch.randn(500,3,4,5)
c = torch.randn(500,1)
#被合并的对象第一维度上的值相等
for x in TensorDataset(b,c):#generator
    print(x)
    break



''''

DataLoader - 用来切割小批量的类

'''


data = TensorDataset(b,c)
from torch.utils.data import DataLoader
for x in DataLoader(data):
    print(x)
    break
bs = 120

dataset = DataLoader(data
          , batch_size = bs
          , shuffle = True #划分小批量之前请随机打乱我们的数据
          , drop_last = True #你要舍弃最后一个batch吗？
          )

for i in dataset:
    print(i[0].shape)

len(dataset) #一共有多少个batch

len(dataset.dataset) #展示里面全部的数据

dataset.dataset[0] #单个样本

dataset.dataset[0][0]

dataset.batch_size #查看现有的batch_size

#TensorDataset - 将特征与标签合并到同一个对象中
#DataLoader - 帮助我们进行小批量的分割









