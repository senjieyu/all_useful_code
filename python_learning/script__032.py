# 反射
#  在程序运行的过程中动态获取程序信息以及动态调用对象的功能

# 为何用反射

# 如何实现反射
class People:
    def __init__(self,name):
        self.name=name

o=People('sad')
print(dir(0))
print(o.__dict__)

#四个内置函数的使用
# hasattr()#判断属性是否存在
print(hasattr(o,'name'))#True

# getattr()   getattr(o,'name',None)
print(getattr(o,'name'))#sad
# setattr()#
setattr(o,'name','jooj')

# delattr()#删除属性
# delattr(o,'name')
o.good=123
print(o.good)




