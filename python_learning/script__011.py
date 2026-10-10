# -*- coding: utf-8 -*-
"""
Created on Wed May 22 10:36:51 2024

@author: senji
"""
import torch
import torch.nn as nn
from torch.nn import functional as F


class Model(nn.Module):
    def __init__(self,in_features=10, out_features=2):
        """
        in_features: 输入该神经网络的特征数目（输入层上的神经元的数目）
        out_features：神经网络的输出的数目（输出层上的神经元的数目）
        """
        super(Model,self).__init__()
        self.linear1 = nn.Linear(in_features,13,bias=True)
        self.linear2 = nn.Linear(13,8,bias=True)
        self.output = nn.Linear(8,out_features,bias=True)
    
    def forward(self,x): #神经网络的向前传播
        z1 = self.linear1(x)
        sigma1 = torch.relu(z1)
        z2 = self.linear2(sigma1)
        sigma2 = torch.sigmoid(z2)
        z3 = self.output(sigma2)
        sigma3 = F.softmax(z3,dim=1)
        return sigma3