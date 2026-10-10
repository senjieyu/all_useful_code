'''
递归锁的特点
可以被第一个抢到锁的 连续的acquire和release
但是只能被第一个抢到这把锁执行上述操作
它的内部有一个计数器 每acquire一次计数加一 每release一次计数减-
只要计数不为0 那么其他人都无法抢到该锁
'''

from threading import  Thread,Lock,RLock
import time

mutexA = mutexB = RLock()

# 类只要加括号多次 产生的肯定是不同的对象
# 如果你想要实现多次加括号等到的是相同的对象   单例模式

class MyThead(Thread):
    def run(self):
        self.func1()
        self.func2()

    def func1(self):
        mutexA.acquire()
        print('%s 抢到A锁'% self.name) # 获取当前线程名
        mutexB.acquire()
        print('%s 抢到B锁'% self.name)
        mutexB.release()
        mutexA.release()

    def func2(self):
        mutexB.acquire()
        print('%s 抢到B锁'% self.name)
        time.sleep(2)
        mutexA.acquire()
        print('%s 抢到A锁'% self.name)
        mutexA.release()
        mutexB.release()

if __name__ == '__main__':
    for i in range(10):
        t = MyThead()
        t.start()











