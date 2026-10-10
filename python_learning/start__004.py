#可以用__file__    当前文件的绝对路径
#print(__file__)
# os.path.dirname(__file__)#拿到__file__的上级路径
#sys.path.append(__file__)
import sys
import os
# print(os.path.dirname(os.path.dirname(__file__)))
PATH1=os.path.dirname(os.path.dirname(__file__))
sys.path.append(PATH1)
from core.src import *

dic[num][1]()  #对输入的num 进行执行





