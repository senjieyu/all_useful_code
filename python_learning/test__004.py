class People:
    school='清华大学'

    def __init__(self,name,sex,age):
        self.name=name
        self.sex=sex
        self.age=age

class Student(People):
    def choose(self):
        print('%s is choosing a course' %self.name)

class Teacher(People):
    def teach(self):
        print('%s is teaching' %self.name)


old=Teacher('asd','asd',45)
print(old.name)
old.teach()
