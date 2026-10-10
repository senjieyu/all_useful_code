# 为类中某个函数加上装饰器@staticmethod后，
# 该函数就变成了非绑定方法，也称为静态方法。
# 该方法不与类或对象绑定，类与对象都可以来调用它，
# 但它就是一个普通函数而已，因而没有自动传值那么一说

import uuid
class MySQL:
    def __init__(self,host,port):
        self.id=self.create_id()
        self.host=host
        self.port=port
    @staticmethod
    def create_id():
        return uuid.uuid1()

conn=MySQL('127.0.0.1',3306)
print(conn.id) #100365f6-8ae0-11e7-a51e-0088653ea1ec

# 类或对象来调用create_id发现都是普通函数，而非绑定到谁的方法
print(MySQL.create_id)
print(conn.create_id)
MySQL.create_id()
conn.create_id()