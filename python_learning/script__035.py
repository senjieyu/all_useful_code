# # 一切皆为对象
# # 什么是元类？
# # 元类就是用来实例化产生类的类
# # 关系：元类------实例化-------》类-----》对象
# class Mymeta(type): #只有继承了type类才能称之为一个元类，否则就是一个普通的自定义类
#     pass
#
# # StanfordTeacher=Mymeta('StanfordTeacher',(object),{...})
# class StanfordTeacher(object,metaclass=Mymeta):
#     school='Stanford'
#
#     def __init__(self,name,age):
#         self.name=name
#         self.age=age
#
#     def say(self):
#         print('%s says welcome to the Stanford to learn Python' %self.name)




class Mymeta(type): #只有继承了type类才能称之为一个元类，否则就是一个普通的自定义类
    def __init__(self,class_name,class_bases,class_dic):
        # print(self) #<class '__main__.StanfordTeacher'>
        # print(class_bases) #(<class 'object'>,)
        # print(class_dic) #{'__module__': '__main__', '__qualname__': 'StanfordTeacher', 'school': 'Stanford', '__init__': <function StanfordTeacher.__init__ at 0x102b95ae8>, 'say': <function StanfordTeacher.say at 0x10621c6a8>}
        super(Mymeta, self).__init__(class_name, class_bases, class_dic)  # 重用父类的功能

        if class_name.islower():
            raise TypeError('类名%s请修改为驼峰体' %class_name)

        if '__doc__' not in class_dic or len(class_dic['__doc__'].strip(' \n')) == 0:
            raise TypeError('类中必须有文档注释，并且文档注释不能为空')

# StanfordTeacher=Mymeta('StanfordTeacher',(object),{...})
class StanfordTeacher(object,metaclass=Mymeta):
    """
    类StanfordTeacher的文档注释
    """
    school='Stanford'

    def __init__(self,name,age):
        self.name=name
        self.age=age

    def say(self):
        print('%s says welcome to the Stanford to learn Python' %self.name)



# 抛出异常 raise