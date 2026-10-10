# 装饰器的内存地址，对于原函数的差别
# index的__name__,__doc__属性不是他自己

# import time
#
# def get_func(func):
#     def wrapper(*args,**kwargs):
#         start=time.time()
#         f1=func(*args,**kwargs)
#         stop=time.time()
#         times=stop-start
#         print(times)
#         return f1
#     #这里修改
#     wrapper().__name__=func().__name__#麻烦，可能会有遗漏
#     wrapper().__doc__=func().__doc__#麻烦，可能会有遗漏
#     return wrapper
#
#
# @get_func
# def index(x,y):
#     ''''这是阿萨大大的文件'''
#     print('输出的是：{}'.format(x))
#     time.sleep(0.1)
#     print('输出的是：{}'.format(y))
#     return 123
#
# index(45,55)


#参数 返回值 属性
#将原函数的属性赋值给wrapper函数
from functools import wraps
import time


def get_func(func):
    #这里将导入的装饰器用一下：@wraps（func）
    @wraps(func)
    def wrapper(*args,**kwargs):
        start=time.time()
        f1=func(*args,**kwargs)
        stop=time.time()
        times=stop-start
        print(times)
        return f1
    return wrapper


@get_func
def index(x,y):
    ''''这是阿萨大大的文件'''
    print('输出的是：{}'.format(x))
    time.sleep(0.1)
    print('输出的是：{}'.format(y))
    return 123

index(45,55)
print(index.__name__)
print(index.__doc__)