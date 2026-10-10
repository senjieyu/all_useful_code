class People:
    def __init__(self,name,weight,height):
        self.name=name
        self.weight=weight
        self.height=height
    @property
    def bmi(self):
        return self.weight / (self.height**2)
obj=People('lili',75,1.85)  #  21.913805697589478
obj.bmi #触发方法bmi的执行，将obj自动传给self，执行后返回值作为本次引用的结果


#property是一个装饰器，用来给绑定对象的方法伪装成一个数据属性

class Foo:
    def __init__(self,val):
        self.__NAME=val #将属性隐藏起来
    @property
    def name(self):
        return self.__NAME
    @name.setter
    def name(self,value):
        if not isinstance(value,str):  #在设定值之前进行类型检查
            raise TypeError('%s must be str' %value)
        self.__NAME=value #通过类型检查后,将值value存放到真实的位置self.__NAME
    @name.deleter
    def name(self):
        raise PermissionError('Can not delete')

f=Foo('lili')
f.name

f.name='LiLi' #触发name.setter装饰器对应的函数name(f,’Egon')
f.name=123 #触发name.setter对应的的函数name(f,123),抛出异常TypeError
del f.name #触发name.deleter对应的函数name(f),抛出异常PermissionError
