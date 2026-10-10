# 在子类派生新方法中如何重用父类的功能
# 方式一：指名道姓调用某一类下的函数===>不依赖继承关系
# class OldboyPenple:
#     def __init__(self,name,age,sex):
#         self.name=name
#         self.age=age
#         self.sex=sex
#
#     def f1(self):
#         print('%s say hello' %self.name)
#
# class Teacher(OldboyPenple):
#     def __init__(self,level,salary):
#         OldboyPenple.__init__(self,name,age,sex)
#         self.level =level
#         self.salary=salary

# 方式二：super（）调用父类提供给自己的方法===>严格依赖继承关系
# 调用super（）会得到一个特殊的对象，该对象会参照属性查找发起者的
# 类的mro，去当前类当前类当前类的父类中找属性
# class OldboyPenple:
#     def __init__(self,name,age,sex):
#         self.name=name
#         self.age=age
#         self.sex=sex
#
#     def f1(self):
#         print('%s say hello' %self.name)
#
# class Teacher(OldboyPenple):
#     def __init__(self,level,salary):
#         # super().__init__(self, name, age, sex)
#         # python3中什么参数都不传，直接加括号就可以
#         super(Teacher,self).__init__(self,name,age,sex)
#         self.level =level
#         self.salary=salary
#


#A没有继承B
class A:
    def test(self):
        print('AAAAAAA')
        super().test()

class B:
    def test(self):
        print('from B')

class C(A,B):
    pass

C.mro() # 在代码层面A并不是B的子类，但从MRO列表来看，属性查找时，就是按照顺序C->A->B->object，B就相当于A的“父类”
# [<class '__main__.C'>, <class '__main__.A'>, <class '__main__.B'>,<class ‘object'>]
obj=C()
obj.test() # 属性查找的发起者是类C的对象obj，所以中途发生的属性查找都是参照C.mro,先找A，A中有super().test()，会参考C.mro去A的父类 <class '__main__.B'>中找test()



#整合
# 在一个类中以另外一个类的对象作为数据属性，称为类的组合。