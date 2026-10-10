import time

def get_func(func):
    def wrapper(*args,**kwargs):
        start=time.time()
        f1=func(*args,**kwargs)
        stop=time.time()
        times=stop-start
        print('运行时间是：',times)
        return f1
    return wrapper

# def p_t(func):
#     def print_s(*args,**kwargs):
#         print('装饰器1尝试{}'.format(args))
#         func(*args,**kwargs)  #返回值改变了
#         print('装饰器1尝试2{}'.format(args))
#     return print_s
def p_t(func):
    def print_s(*args,**kwargs):
        print('装饰器1尝试{}'.format(args))
        f1=func(*args,**kwargs)  #返回值改变了
        print('装饰器1尝试2{}'.format(args))
        return f1
    return print_s

@p_t
@get_func  #写装饰器的名字
def index(x,y):
    print('输出的是：{}'.format(x))
    time.sleep(0.1)
    print('输出的是：{}'.format(y))
    return 123

@p_t
@get_func
def good_boy():
    print("asdd")
    print('zheshire')
    return 456

f1=index(6,8)
print(f1)


