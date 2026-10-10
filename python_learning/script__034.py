# class People:
#     def __init__(self,name,age):
#         self.name=name
#         self.age=age
#     def __str__(self):
#         return '<Name:%s Age:%s>' %(self.name,self.age) #返回类型必须是字符串
# #__str__：在打印对象时自动触发，然后将返回值（必须是字符类型）当做本次打印的结果输出
# p=People('lili',18)
# print(p) #触发p.__str__()，拿到返回值后进行打印


# __del__会在对象被清理时自动触发。
# 由于Python自带的垃圾回收机制会自动清
# 理Python程序的资源，所以当一个对象只占用
# 应用程序级资源时，完全没必要为对象定制__del__方法，
# 但在产生一个对象的同时涉及到申请系统资源
# （比如系统打开的文件、网络连接等）的情况下，
# 关于系统资源的回收，Python的垃圾回收机制便派不上用场了，
# 需要我们为对象定制该方法，用来在对象被删除时自动触发回收系统资源的操作

# class MySQL:
#     def __init__(self,ip,port):
#         self.conn= connect(ip,port) # 伪代码，发起网络连接，需要占用系统资源
#     def __del__(self):
#         self.conn.close() # 关闭网络连接，回收系统资源
#
# obj=MySQL('127.0.0.1',3306) # 在对象obj被删除时，自动触发obj.__del__()






