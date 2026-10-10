# class ATM:
#     def __card(self): #插卡
#         print('插卡')
#     def __auth(self): #身份认证
#         print('用户认证')
#     def __input(self): #输入金额
#         print('输入取款金额')
#     def __print_bill(self): #打印小票
#         print('打印账单')
#     def __take_money(self): #取钱
#         print('取款')
#     def withdraw(self): #取款功能
#         self.__card()
#         self.__auth()
#         self.__input()
#         self.__print_bill()
#         self.__take_money()
# obj = ATM()
# obj.withdraw()

# class People:
#     def __init__(self,name,weight,height):
#         self.name=name
#         self.weight=weight
#         self.height=height
#     @property
#     def bmi(self):
#         return self.weight / (self.height**2)
#
# obj=People('lili',75,1.85)
# print(obj.bmi) #触发方法bmi的执行，将obj自动传给self，执行后返回值作为本次引用的结果




# 派生
# class People:
#     school='清华大学'
#
#     def __init__(self,name,sex,age):
#         self.name=name
#         self.sex=sex
#         self.age=age
#
# # class Teacher(People):
# #     def __init__(self,name,sex,age,title): # 派生
# #         self.name=name
# #         self.sex=sex
# #         self.age=age
# #         self.title=title
# #     def teach(self):
# #         print('%s is teaching' %self.name)
# #
# # obj=Teacher('lili','female',28,'高级讲师') #只会找自己类中的__init__，并不会自动调用父类的
#
# # 方法1
# # class Teacher(People):
# #     def __init__(self,name,sex,age,title):
# #         People.__init__(self,name,age,sex) #调用的是函数,因而需要传入self
# #         self.title=title
# #     def teach(self):
# #         print('%s is teaching' %self.name)
# #
# # 方法2
# # class Teacher(People):
# #     def __init__(self,name,sex,age,title):
# #         super().__init__(name,age,sex) #调用的是绑定方法，自动传入self
# #         self.title=title
# #     def teach(self):
# #         print('%s is teaching' %self.name)
# 这两种方式的区别是：方式一是跟继承没有关系的，
# 而方式二的super()是依赖于继承的，并且即使没有直接继承关系，
# super()仍然会按照MRO继续往后查找





# # 类的组合
# class Course:
#     def __init__(self,name,period,price):
#         self.name=name
#         self.period=period
#         self.price=price
#     def tell_info(self):
#         print('<%s %s %s>' %(self.name,self.period,self.price))
#
# class Date:
#     def __init__(self,year,mon,day):
#         self.year=year
#         self.mon=mon
#         self.day=day
#     def tell_birth(self):
#        print('<%s-%s-%s>' %(self.year,self.mon,self.day))
#
# class People:
#     school='清华大学'
#     def __init__(self,name,sex,age):
#         self.name=name
#         self.sex=sex
#         self.age=age
#
# #Teacher类基于继承来重用People的代码，基于组合来重用Date类和Course类的代码
# class Teacher(People): #老师是人
#     def __init__(self,name,sex,age,title,year,mon,day):
#         super().__init__(name,age,sex)
#         self.birth=Date(year,mon,day) #老师有生日
#         self.courses=[] #老师有课程，可以在实例化后，往该列表中添加Course类的对象
#     def teach(self):
#         print('%s is teaching' %self.name)
#
#
# python=Course('python','3mons',3000.0)
# linux=Course('linux','5mons',5000.0)
# teacher1=Teacher('lili','female',28,'博士生导师',1990,3,23)
#
# # teacher1有两门课程
# teacher1.courses.append(python)
# teacher1.courses.append(linux)
#
# # 重用Date类的功能
# teacher1.birth.tell_birth()
#
# # 重用Course类的功能
# for obj in teacher1.courses:
#     obj.tell_info()



#
# # 多态性
# class Animal: #同一类事物:动物
#     def talk(self):
#         pass
# class Cat(Animal): #动物的形态之一:猫
#     def talk(self):
#         print('喵喵喵')
# class Dog(Animal): #动物的形态之二:狗
#     def talk(self):
#         print('汪汪汪')
# class Pig(Animal): #动物的形态之三:猪
#     def talk(self):
#         print('哼哼哼')
#
# #实例化得到三个对象
# cat=Cat()
# dog=Dog()
# pig=Pig()
# # 多态性指的是可以在不用考虑对象具体类型的情况下而直接使用对象，
# # 这就需要在设计时，把对象的使用方法统一成一种：例如cat、dog、pig都是动物,
# # 但凡是动物肯定有talk方法，
# # 于是我们可以不用考虑它们三者的具体是什么类型的动物,而直接使用
# # 于是我们可以不用考虑它们三者的具体是什么类型的动物,而直接使用
# # 于是我们可以不用考虑它们三者的具体是什么类型的动物,而直接使用
# # 于是我们可以不用考虑它们三者的具体是什么类型的动物,而直接使用
# # 于是我们可以不用考虑它们三者的具体是什么类型的动物,而直接使用
# # 于是我们可以不用考虑它们三者的具体是什么类型的动物,而直接使用
# cat.talk()
# dog.talk()
# pig.talk()
# # 更进一步,我们可以定义一个统一的接口来使用
# def Talk(animal):
#     animal.talk()
#
# Talk(cat)
# Talk(dog)
# Talk(pig)
#
#
# # Python中一切皆对象，本身就支持多态性
# # 我们可以在不考虑三者类型的情况下直接使用统计三个对象的长度
# 's'.__len__()
# [45,55].__len__()
# ('ad','t').__len__()
#
# # Python内置了一个统一的接口
# len('s')
# len([45,55])
# len(('ad','t'))
#
# # 综上我们得知，多态性的本质在于不同的类中定义有相同的方法名，
# # 这样我们就可以不考虑类而统一用一种方式去使用对象，
# # 可以通过在父类引入抽象类的概念来硬性限制子类必须有某些方法名



# # 综上我们得知，多态性的本质在于不同的类中定义有相同的方法名，
# # 这样我们就可以不考虑类而统一用一种方式去使用对象，
# # 可以通过在父类引入抽象类的概念来硬性限制子类必须有某些方法名
# # 鸭子类型
# # 二者看起来都像文件,因而就可以当文件一样去用，然而它们并没有直接的关系
# class Txt: # Txt类有两个与文件类型同名的方法，即read和write
#     def read(self):
#         pass
#     def write(self):
#         pass
#
# class Disk: # Disk类也有两个与文件类型同名的方法：read和write
#     def read(self):
#         pass
#     def write(self):
#         pass
#









