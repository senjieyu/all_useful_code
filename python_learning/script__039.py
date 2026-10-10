'''信号量在不同的阶段可能对应不同的技术点
在并发编程中信号量指的是锁!!!
如果我们将互斥锁比喻成一个厕所的话
那么信号量就相当于多个厕所

           多把锁
'''

from threading import Thread,Semaphore
import time
import random

'''
利用random模块实现打印随机验证码
'''
sm = Semaphore(5)   # 括号内写几  就代表几个坑位

def task(name):
    sm.acquire()
    print('%s 正在蹲坑'%name)
    time.sleep(random.randint(1,5))
    sm.release()

if __name__ == '__main__':
    for i in range(20):
        t = Thread(target=task,args = ('伞兵%s号'%i, ))
        t.start()






