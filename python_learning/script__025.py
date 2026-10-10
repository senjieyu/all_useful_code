from multiprocessing import Process
from threading import Thread
import time


class MyThead(Thread):
    def __init__(self,name):
        super().__init__()
        self.name = name
    def run (self):
        print('%s is running'%self.name)
        time.sleep(1)
        print('1111111')

if __name__ == '__main__':



    # 开启线程不需要再main下面执行代码 直接书写就可以
    # 但是我们还是习惯性的将启动命令写在main下面
    t = MyThead('sad')
    t.start()
    print('主')
