# 封装是面向对象三大特性最核心的特性
# 分装《<--->》整合
#


#  二、将封装的属性进行隐藏操作
# 如何隐藏： 在属性名前加__前缀，就会实现对外隐藏的效果
# 该隐藏需要注意的问题
# 1.在类外部无法直接访问双下滑线开头的属性，但知道了类名和属性名就可
# 以拼出名字：_类名__属性，然后就可以访问了，如Foo._A__N，所以说
# 这种操作并没有严格意义上地限制外部访问，仅仅只是一种语法意义上的变形。
# class Foo:
#     __x=1   # _Foo__x     只是一种变形操作
#
#     def __f1 (self):
#         print('test')

# 2.这种隐藏对外不对内,因为__开头的属性会在检查类内
# 代码语法时统一发生变形

# 3.这种变形操作只在检查类语法的时候发生一次，之后定义的__开头的属性都不会发生变形

# class Foo:
#     __x=1   # _Foo__x     只是一种变形操作
#
#     def __f1 (self):
#         print('test')
#
#     def f2 (self):
#         print(self.__x)
#         print(self.__f1())
# obj=Foo()
# obj.f2()
# Foo.__y=3
# print(Foo.__y) # __y  不变形


#为何要隐藏？
#隐藏数据属性
#隐藏函数/方法属性

class People:
    def __init__(self,name):
        self.__name=name

    def get_name(self):
        #通过该 接口 可以间接的访问到名字属性
        print(self.name)


odj=People('sad')
odj.get_name()

# 将数据隐藏起来就限制了类外部对数据的直接操作，
# 然后类内应该提供相应的接口来允许类外部间接地操
# 作数据，接口之上可以附加额外的逻辑来对数据的操
# 作进行严格地控制


#
# 目的的是为了  隔离复杂度  ，例如ATM程序的取款功能,该功
# 能有很多其他功能组成，比如插卡、身份认证、输入金额、
# 打印小票、取钱等，而对使用者来说,只需要开发取款这个
# 功能接口即可,其余功能我们都可以隐藏起来
class ATM:
    def __card(self): #插卡
        print('插卡')
    def __auth(self): #身份认证
        print('用户认证')
    def __input(self): #输入金额
        print('输入取款金额')
    def __print_bill(self): #打印小票
        print('打印账单')
    def __take_money(self): #取钱
        print('取款')
    def withdraw(self): #取款功能
        self.__card()
        self.__auth()
        self.__input()
        self.__print_bill()
        self.__take_money()
obj=ATM()
obj.withdraw()

