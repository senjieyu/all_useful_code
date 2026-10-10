from multiprocessing import Process
from threading import Thread
import time


def task(name):
    print('%s is running'%name,end = '\n')
    time.sleep(1)
    print('%s is over'%name)


# 开启线程不需要再main下面执行代码 直接书写就可以
# 但是我们还是习惯性的将启动命令写在main下面
t = Thread(target=task, args=('二狗',))    # 线程不需要申请内存空间
# p = Process(target=task, args=('erdan',))  进程会先申请内存空间
# 而主进程代码会继续以高速往下运行，所以会先打印 主
t.start()
print('主')




























