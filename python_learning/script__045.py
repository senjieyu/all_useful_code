import time

def get_func(func):
    def wrapper(*args,**kwargs):
        start=time.time()
        f1=func(*args,**kwargs)
        stop=time.time()
        times=stop-start
        print(times)
        return f1
    return wrapper

@get_func  #写装饰器的名字
def index(x,y):
    print('输出的是：{}'.format(x))
    time.sleep(0.1)
    print('输出的是：{}'.format(y))
    return 123

@get_func
def good_boy():
    print("asdd")
    print('zheshire')


f1=index(45,46)
print(f1)
good_boy()


#叠加多个装饰器
# @sada
# @asdas
# @asddas
# def index
#运行顺序：自上而下，装饰顺序自下而上

