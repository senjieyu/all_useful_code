# # 进程间通信
from multiprocessing import Queue, Process

# '''
# 1.主进程跟子进程借助队列通信
#
# '''
# def producer(q):
#     print('qweasdzxc')
#     q.put('我打算')
#
# if __name__ == '__main__':
#     q = Queue()
#     p = Process(target = producer , args=(q,))
#     p.start()
#     print(q.get())

'''
1.子进程跟子进程借助
'''
import time

def producer(q):
    q.put('我打算1')


def consumer(q):
    print(q.get())


if __name__ == '__main__':
    q = Queue()
    p = Process(target=producer, args=(q,))
    p1 = Process(target=consumer, args=(q,))
    p1.start()
    time.sleep(1)
    p.start()
