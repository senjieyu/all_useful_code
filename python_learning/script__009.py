# -*- coding: utf-8 -*-
"""
Created on Thu May 23 09:13:36 2024

@author: senji
"""

# #带大家认识fashion-MINST数据集
# 本节课我们讲解了神经网络使用小批量随机梯度下降进行迭代的流程，现在我们要整合本节课中所有的
# 代码实现一个完整的训练流程。首先要梳理一下整个流程：
# 1）设置步长 ，动量值 ，迭代次数 ，batch_size等信息，（如果需要）设置初始权重
# ，
# 2）导入数据，将数据切分成batches
# 3）定义神经网络架构
# 4）定义损失函数 ，如果需要的话，将损失函数调整成凸函数，以便求解最小值
# 5）定义所使用的优化算法
# 6）开始在epoches和batch上循环，执行优化算法：
# 6.1）调整数据结构，确定数据能够在神经网络、损失函数和优化算法中顺利运行
# 6.2）完成向前传播，计算初始损失
# 6.3）利用反向传播，在损失函数 上对每一个 求偏导数
# 6.4）迭代当前权重
# 6.5）清空本轮梯度
# 6.6）完成模型进度与效果监控
# 7）输出结果

import torch
from torch import nn
from torch import optim
from torch.nn import functional as F
from torch.utils.data import DataLoader, TensorDataset
import torchvision
import torchvision.transforms as transforms 
import torchvision
import torchvision.transforms as transforms #处理数据模块
#dataloader、TensorDataset - 对数据的结构、归纳方式进行变换
#torchvision.transforms - 对数据集的数字本身进行修改

mnist = torchvision.datasets.FashionMNIST(root = "D:/pytorch_learning/MINST-FASHION" #你的计算机上的某个目录
                                          ,download = False
                                          , train = True
                                          , transform = transforms.ToTensor()
                                         ) #实例化数据

#TensorDataset 不需要
#小批量的划分

#确定数据、确定超参数

lr = 0.15
gamma = 0
epochs = 10
bs = 128

batchdata = DataLoader(mnist
                       ,batch_size = bs
                       ,shuffle = True)

for x,y in batchdata:
    print(x.shape)
    print(y.shape)
    break
# torch.Size([128, 1, 28, 28])
# torch.Size([128])


input_ = mnist.data[0].numel() #请问这个张量中总共有多少元素呢？
output_ = len(mnist.targets.unique())



class Model(nn.Module):
    def __init__(self, in_features=10, out_features=2):
        super().__init__()
        self.linear1 = nn.Linear(in_features,1280,bias=False)
        self.output = nn.Linear(1280,out_features, bias=False)
    
    def forward(self, x):
        x = x.view(-1,28*28)
        
        # 它类似于 NumPy 中的 reshape 方法, view(-1,) #需要对数据结构进行一个改变，-1作为占位符，表示请pytorch帮助我们自动计算-1这个位置的维度应是多少
        '''
        x = torch.randn(30,40)
        x = x.view(-1,20)
        x.shape
        torch.Size([60, 20])
        '''
        
        sigma1 = torch.relu(self.linear1(x))
        sigma2 = F.log_softmax(self.output(sigma1),dim=1) # 多分类中softmax
        return sigma2


def fit(net,bacthdata,lr=0.01, epochs=5, gamma = 0):
    criterion = nn.NLLLoss() #定义损失函数
    opt = optim.SGD(net.parameters(),lr=lr,momentum = gamma) #小批量梯度下降
    correct = 0 # 循环开始之前，预测正确的值为0
    samples = 0 # 循环开始之前，模型一个样本都没有见过
    for epoch in range(epochs): # 全数据被训练几次
        for batch_idx,(x,y) in enumerate(batchdata):   # enumerate   对迭代器对象进行标号
            '''list = ["A","B","C"]
            for x in enumerate(list):
                print(x)
            (0, 'A')
            (1, 'B')
            (2, 'C')
            '''
            y = y.view(x.shape[0]) #降维  保证y为1维
            sigma = net.forward(x) #正向传播
            loss = criterion(sigma,y)  # 损失函数 
            loss.backward()    #反向传播
            opt.step()    #更新权重
            opt.zero_grad()    #梯度清零
            
            #求解准确率，全部判断正确的样本数量/已经看过的总样本量
            yhat = torch.max(sigma, 1)[1] # torch.max函数结果中的索引为1的部分     用于获取张量中的最大值。它有多种用法，具体取决于传递的参数。常见的用法包括查找整个张量中的最大值，以及沿指定维度查找最大值和相应的索引。
            correct += torch.sum(yhat == y)
            samples += x.shape[0]
            #每训练一个batch的数据，模型见过的数据就会增加x.shape[0]
            
            if (batch_idx + 1) % 125 == 0 or batch_idx == len(batchdata)-1: #每125个batch我就打印一次
                print("Epoch{}:[{}/{}({:.0f}%)],Loss:{:.6f},Accuracy:{:.3f}".format(
                    epoch+1
                    ,samples
                    ,epochs*len(batchdata.dataset)
                    ,100*samples/(epochs*len(batchdata.dataset))
                    ,loss.data.item()
                    ,float(100*correct/samples)
                     ))
            # 分子代表：已经查看过的数据有多少
            # 分母代表：在现有的epochs设置，模型一共需要查看多少数据


#训练与评估
torch.manual_seed(420)
net = Model(in_features=input_, out_features=output_)
fit(net,batchdata,lr=lr,epochs=epochs,gamma=gamma)




















