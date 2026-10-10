# 同一个进程下多个线程数据是共享的
from threading import Thread
import time
money = 100

def task():
    global  money
    money = 666

if __name__ == '__main__':
    t = Thread(target=task)
    t.start()
    print(money)







