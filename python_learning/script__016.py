'''
#僵尸进程
死了但没有死透
当你开设了子进程之后 该进程死后不会立刻释放占用的进程号
因为我要让父进程能够查到它开设的子进程的一些基本信息 占用的pid号 运行时间
所有的进程都会步入僵尸进程
    父进程不死并且在无限制的创建子进程并且子进程也不结束
    回收进程占用的pid号
    父进程等待子进程运行结束
    父进程调用join方法
'''

import time
from multiprocessing import Process

def run():
    print('hello world')
    time.sleep(3)
    print('get_out')

if __name__ == '__main__':
    p = Process(target=run)
    p.start()
    print('zhu')

"""
孤儿进程
子进程存活，父进程意外死亡
操作系统会开设一个’儿童福利院‘专门管理孤儿进程回收相关资源
"""