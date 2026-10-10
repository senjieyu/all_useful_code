# -*- coding: utf-8 -*-
"""
Created on Wed May 22 16:52:02 2024

@author: senji
"""
##导入库
#确定数据、超参数的确定（lr，gamma）
#定义伸进网络的架构类Model，类Model需要输入的参数
##实例化神经网络的类 - 让神经网络准备好进行正向传播
#定义损失函数
#定义优化算法

import torch
import torch.nn as nn
import torch.optim as optim
from torch.nn import functional as F
#确定数据
torch.manual_seed(420)
X = torch.rand((500,20),dtype=torch.float32) * 100
y = torch.randint(low=0,high=3,size=(500,),dtype=torch.float32)

lr = 0.1
gamma = 0.9


#定义神经网路的架构
"""
注意：这是一个三分类的神经网络，因此我们需要调用的损失函数多分类交叉熵函数CEL
CEL类已经内置了softmax功能，因此我们需要修改一下网络架构，删除forward函数中输出层上的softmax函数，并将最终的输出修改为zhat
"""
class Model(nn.Module):
    def __init__(self,in_features=10,out_features=2):
        super().__init__() #super(请查找这个类的父类，请使用找到的父类替换现在的类)
        self.linear1 = nn.Linear(in_features,13,bias=True) #输入层不用写，这里是隐藏层的第一层
        self.linear2 = nn.Linear(13,8,bias=True)
        self.output = nn.Linear(8,out_features,bias=True)
        
    def forward(self, x):
        sigma1 = torch.relu(self.linear1(x))
        sigma2 = torch.sigmoid(self.linear2(sigma1))
        zhat = self.output(sigma2)
        return zhat


input_ = X.shape[1] #特征的数目
output_ = len(y.unique()) #分类的数目

#实例化神经网络类
torch.manual_seed(420)
net = Model(in_features=input_, out_features=output_)

net.linear1.weight
net.parameters()
for x in net.parameters():
    print(x)

#定义损失函数
criterion = nn.CrossEntropyLoss()

opt = optim.SGD(net.parameters()#需要进行迭代的权重
                ,lr = lr
                ,momentum = gamma
               )

#向前传播
#本轮向前传播的损失函数值
#反向传播 - 得到了梯度
#更新权重（和动量）
#清空梯度 - 清除原来计算出来的，基于上一个点的坐标计算的梯度


zhat = net.forward(X) #最后一个线性层的输出结果，向前传播
loss = criterion(zhat, y.reshape(500).long()) #计算损失函数
loss.backward()
opt.step() #步子，走一步，更新权重w，更新动量v
opt.zero_grad()

print(loss)
print(net.linear1.weight.data[0][:10])
















