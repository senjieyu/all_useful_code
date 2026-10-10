# import time
# def index(x,y):
#     print('输出的是：{}'.format(x))
#     time.sleep(0.1)
#     print('输出的是：{}'.format(y))
#     return 123
# def get_func(func):
#     def wrapper(*args,**kwargs):
#         start=time.time()
#         func(*args,**kwargs)
#         stop=time.time()
#         times=stop-start
#         print(times)
#         return func(*args,**kwargs)   #这里又调用了一次
#     return wrapper
# index=get_func(index)
# f1=index(48,18)
# print(f1)

import time
def index(x,y):
    print('输出的是：{}'.format(x))
    time.sleep(0.1)
    print('输出的是：{}'.format(y))
    return 123
def get_func(func):
    def wrapper(*args,**kwargs):
        start=time.time()
        f1=func(*args,**kwargs)
        stop=time.time()
        times=stop-start
        print(times)
        return f1
    return wrapper
index=get_func(index)
f1=index(5575,1757)
print(f1)