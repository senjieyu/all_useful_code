# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:12:14 2024

@author: senji
"""

import torch
import torch.nn as nn
from torch.nn import CrossEntropyLoss as CEL
from torch.nn import functional as F

torch.manual_seed(420)
X = torch.rand((500,20),dtype=torch.float32) * 100
y = torch.randint(low=0,high=3,size=(500,),dtype=torch.float32)

input_ = X.shape[1]
output_ = len(y.unique())

class Model(nn.Module):
    def __init__(self,in_features=40,out_features=2):
        super().__init__()
        self.linear1 = nn.Linear(in_features,13,bias=False)
        self.linear2 = nn.Linear(13,8,bias=False)
        self.output = nn.Linear(8,out_features,bias=True)
    
    def forward(self, x):
        sigma1 = torch.relu(self.linear1(x))
        sigma2 = torch.sigmoid(self.linear2(sigma1))
        zhat = self.output(sigma2)
        return zhat

torch.manual_seed(420)
net = Model(in_features=input_, out_features=output_)

zhat = net.forward(X)

criterion = CEL()

loss = criterion(zhat, y.long())


loss.backward(retain_graph=True)


print(net.linear1.weight.grad)
















