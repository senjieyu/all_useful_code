class Vehicle:  # 交通工具
    pass


class FlyableMixin:
    def fly(self):
        '''
        飞行功能相应的代码
        '''
        print("I am flying")


class CivilAircraft(FlyableMixin, Vehicle):  # 民航飞机
    pass


class Helicopter(FlyableMixin, Vehicle):  # 直升飞机
    pass


class Car(Vehicle):  # 汽车
    pass

# ps: 采用某种规范（如命名规范）来解决具体的问题是python惯用的套路