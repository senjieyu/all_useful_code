# 1.什么是继承
#继承是一种创建新类的方式，在Python中，新
# # 建的类可以继承一个或多个父类，新建的类可
# # 称为子类或派生类，父类又可称为基类或超类
# # python支持多继承
# # 在python中，新建的类可以继承一个或多个父类
#
# class Parent1:
#     x=111
#
# class Parent2:
#     pass
#
# class Sub1(Parent1):#单继承
#     pass
#
# class Sub2(Parent1,Parent2):#多继承
#     pass
#
# print(Sub1.__bases__)#查看父类 (<class '__main__.Parent1'>,)
# print(Sub2.__bases__)#(<class '__main__.Parent1'>, <class '__main__.Parent2'>)
# # ps：在python2中有经典类与新式类之分
# # 新式类：继承了object类的子类，以及该子类的子类的子类。。
# # 经典：没有继承object类的子类，以及该子类的子类子子类。。
# print(Parent1.__bases__)#(<class 'object'>,)
# print(Parent2.__bases__)#(<class 'object'>,)
# # 在python中没有继承任何类，那么会默认继
# # 承object类，所以python3中所有的类都是新式类  object
# # object类是python的内置类，包含一些比较重要的功能
#
#
# print(Sub1.x)  #  111
#
#
#
# # python的多继承
# # 优点：子类可以遗传多个父类的属性，最大限度地重用代码
# # 缺点：
# #     违背人的思维习惯：继承表达的是一种什么‘是’什么的关系
# #     代码的可读性变差
# #     不建议使用多继承，有可能会引发菱形问题，扩展性变差，应该用Mixins
#
# 2.为何要用继承：
# 子类会遗传父类的属性，可以减少类与类之间代码的冗余
# 用来解决代码冗余的问题
# 派生,自己有则以自己为准


# 3.如何实现继承
# class Student:
#     school='asd'
#
#     def __int__(self,name,age,sex):
#         self.name=name
#         self.age=age
#         self.sex=sex
#
#     def choose_course (self):
#         print('%s 行啊额恶化' %self.name)
#
#
# class Teacher:
#     school='asdf'
#     def __int__(self, name, age, sex,salary,level):
#         self.name = name
#         self.age = age
#         self.sex = sex
#         self.salary = salary
#         self.level = level
#     def score(self):
#         print('%s老师正在教学' %self.name)


# 类与类之间存在冗余
# class Oldboy:
#     school = 'sadsasf'
#     def __init__(self, name, age, sex):
#         self.name = name
#         self.age = age
#         self.sex = sex
#
# class Student(Oldboy):
#     school = 'asd45455'
#     def choose_course (self):
#         print('%s 行啊额恶化' %self.name)
#
# class Teacher(Oldboy):
#     def __init__(self, name, age, sex,salary,level):
#         Oldboy.__init__(self, name, age, sex)
#         self.salary = salary
#         self.level = level
#     def score(self):
#         print('%s老师正在教学' %self.name)
#
# test = Teacher('asd','1','3','4','5')
# print(test.name)
# print(test.school)



# class Foo:
#     def f1(self):
#         print('Foo.f1')
#     def f2(self):
#         print('Foo.f2')
#         self.f1()
#
# class Bar(Foo):
#     def f1(self):
#         print('Bar.f1')
#
# b=Bar()
# b.f2()
# # Foo.f2
# # Bar.f1   先从自己找



class Foo:
    def __f1(self): # 变形为_Foo__fa
        print('Foo.f1')
    def f2(self):
        print('Foo.f2')
        self.__f1() # 变形为self._Foo__fa,因而只会调用自己所在的类中的方法

class Bar(Foo):
    def __f1(self): # 变形为_Bar__f1
        print('Bar.f1')

b=Bar()
b.f2() #在父类中找到f2方法，进而调用b._Foo__f1()方法，同样是在父类中找到该方法
# Foo.f2
# Foo.f1