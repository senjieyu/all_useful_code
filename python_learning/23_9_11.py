import time
def index(x,y):
    print('输出的是：{}'.format(x))
    time.sleep(3)
    print('输出的是：{}'.format(y))
    print(f'asdasdas{x}')

def wrapper(*args,**kwargs):
    start=time.time()
    index(*args,**kwargs)
    stop=time.time()
    times=stop-start
    print(times)



wrapper(11,22)