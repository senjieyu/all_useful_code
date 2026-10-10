# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:38:54 2024

@author: senji
"""
#导入库、数据、定义神经网络类，完成正向传播

#继承nn.Module类完成正向传播
import torch
import torch.nn as nn
from torch.nn import functional as F

#确定数据
torch.manual_seed(420)
X = torch.rand((500,20),dtype=torch.float32) * 100
y = torch.randint(low=0,high=3,size=(500,),dtype=torch.float32)

#定义神经网路的架构
"""
注意：这是一个三分类的神经网络，因此我们需要调用的损失函数多分类交叉熵函数CEL
CEL类已经内置了softmax功能，因此我们需要修改一下网络架构，删除forward函数中输出层上的softmax函数，并将最终的输出修改为zhat
"""
class Model(nn.Module):
    def __init__(self,in_features=10,out_features=2):
        super(Model,self).__init__() #super(请查找这个类的父类，请使用找到的父类替换现在的类)
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

#前向传播
zhat = net.forward(X)

#定义损失函数
criterion = nn.CrossEntropyLoss()
#对打包好的CorssEnrtopyLoss而言，只需要输入zhat
loss = criterion(zhat,y.long())

loss

net.linear1.weight.grad #不会返回任何值

#反向传播，backward是任意损失函数类都可以调用的方法，对任意损失函数，backward都会求解其中全部w的梯度
loss.backward()

net.linear1.weight.grad #返回相应的梯度

#与可以重复进行的正向传播不同，一次正向传播后，反向传播只能进行一次
#如果希望能够重复进行反向传播，可以在进行第一次反向传播的时候加上参数retain_graph
loss.backward(retain_graph=True)



lr = 0.1
gamma = 0.9#learning_rate, 0.001,0.01,0.05
w = net.linear1.weight.data
dw = net.linear1.weight.grad
w -= lr*dw

#momentum
# v(t) = gamma * v(t-1)  - lr * dw
# w(t+1) = w(t) + v(t)


# 其中，对上一步的梯度向量加上的权重被称为动量参数（也叫做衰减力度，通常使用 进行表示），对这
# 一点的梯度向量加上的权重就是步长（依然表示为 ），真实移动的向量为 ，被称为”动量
# “（Momentum）


#t = 1，走第一步，进行首次迭代的时候，需要一个v0


v = torch.zeros(dw.shape[0], dw.shape[1])

v = gamma*v - lr*dw
w += v
































