class Studen:
    def __init__(self,name,age,stage,couse):
        self.name=name
        self.age=age
        self.stage=stage
        self.couse = couse
    def check_stu(self):
        print('%s学生今年%d岁%d年级'%(self.name,self.age,self.stage))


class Teature:
    def __init__(self, name, age, couse):
        self.name = name
        self.age = age
        self.couse=couse


class ClassCourse:
    def __init__(self,shool,time,money):
        self.shool = shool
        self.time = time
        self.money=money

class School:
    def __init__(self,shool,time,money):
        self.shool = shool
        self.time = time
        self.money=money

stu1=Studen('wang',15,9)
stu1.check_stu()

