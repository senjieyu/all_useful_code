import sys
sys.path.append(r"C:\Users\senji\Desktop\python_learning\learning\23.9.16包的尝试\foo")
#可以用__file__    当前文件的绝对路径
#print(__file__)
#sys.path.append(__file__)
# import os
# os.path.dirname(__file__)#拿到__file__的上级路径
print(sys.path)
from foo import *
f6()
f2()



