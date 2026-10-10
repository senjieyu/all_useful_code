#绑定方法：特殊之处在于将调用者本身当做第一个参数自动传入
# 1.绑定给对象的方法：调用者是对象，自动传入的是对象
# # 2.绑定给类的方法：调用者是类，自动传入的是类
# # 配置文件settings.py的内容
# HOST='127.0.0.1'
# PORT=3306

# @classmethod装饰器的应用

# 类方法的应用
# import settings
# class MySQL:
#     def __init__(self,host,port):
#         self.host=host
#         self.port=port
#     @classmethod   #将下面的函数装饰成绑定给类的方法
#     def from_conf(cls): # 从配置文件中读取配置进行初始化
#         return cls(settings.HOST,settings.PORT)
#
# >>> MySQL.from_conf # 绑定到类的方法
# <bound method MySQL.from_conf of <class ‘__main__.MySQL'>>
# >>> conn=MySQL.from_conf() # 调用类方法，自动将类MySQL当作第一个参数传给cls




#非绑定方法




