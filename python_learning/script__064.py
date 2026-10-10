# 多态指的是一类事物有多种形态， 比如动物有多种形态：猫、狗、猪

class Animal: #同一类事物:动物
    def talk(self):
        print('会叫')
class Cat(Animal): #动物的形态之一:猫
    def talk(self):
        super().talk()
        print('喵喵喵')
class Dog(Animal): #动物的形态之二:狗
    def talk(self):
        super().talk()
        print('汪汪汪')
class Pig(Animal): #动物的形态之三:猪
    def talk(self):
        super().talk()
        print('哼哼哼')

#为何要有多态，多态会带来什么样的特性，多态性

obj1=Cat()
obj2=Dog()
obj3=Pig()

#定义统一接口
def animal_say(animal):
    animal.talk()

animal_say(obj3)
animal_say(obj1)
animal_say(obj2)





