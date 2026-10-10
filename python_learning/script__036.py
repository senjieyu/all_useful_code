# 先回顾之前TCP服务端实现并发的效果是怎么玩的
# 每来一个人就开设一个进程或者线程去处理
# 无论是开设进程也好还是开设线程也好 是不是都需要消耗资源
# 只不过开设线程的消耗比开设进程的稍微小一点而已
# 我们是不可能做到无限制的开设进程和线程的 因为计算机硬件的资源更不上!!!
# 硬件的开发速度远远赶不上软件

# 我们的宗旨应该是保证计算机硬件能够正常工作的限度下最大限度的利用它

# 池的概念
"""
什么是池？
池是用来保证计算机硬件安全的情况下最大限度的利用计算机
它降低了程序的运行效率但是保证了计算机硬件的安全 从而让你写的程序能够正常运行
"""

from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time,os



# pool = ThreadPoolExecutor(5)  # 池子里固定只有5个线程
# 括号内可以传数字 不传的话默认会开设当前计算机cpu个数五倍的线程

pool = ProcessPoolExecutor()  # 池子里固定只有5个线程
# 括号内可以传数字 不传的话默认会开设当前计算机cpu个数的进程



'''
池子造出来之后 里面会固定存在五个线程
这个五个线程不会出现重复创建和销毁的过程

池子造出来之后 里面会固定存在几个进程
这个几个进程不会出现重复创建和销毁的过程

池子的使用非常的简单
你只需要将需要做的任务往池子中提交即可 自动会有人来服务你
'''

def task(n):
    print(n,os.getpid())
    time.sleep(2)
    return n**n

def call_back(n):
    print('call_back>>>:',n.result())

'''
任务的提交方式
同步:提交任务之后原地等待任务的返回结果 期间不做任何事
异步:提交任务之后不等待任务的返回结果 直接继续往下执行
     返回结果如何获取？？？     
     
     异步提交任务的返回结果 应该通过  回调函数  来获取
     一旦该任务有结果立刻触发爆炸
'''
if __name__ == '__main__':
    #  进程必须在main下面，不然会报错

    # pool.submit(task,1) # 朝池子中提交任务   异步提交
    # print('主')
    t_list = []
    for i in range(20):
        # res = pool.submit(task,i)
        # 异步回调机制
        res = pool.submit(task,i).add_done_callback(call_back)

        # print(res.result()) # result方法   同步提交  join
        t_list.append(res)
    # # 等待线程池中所有的任务执行完毕之后再继续往下执行
    # pool.shutdown()  # 关闭线程池  等待线程池中所有的任务运行完毕
    #
    #
    # for t in t_list:
    #     print('>>>:',t.result())   # 肯定是有序的
'''
.result方法
程序由并发变成了串行
任务的为什么打印的是None
res.result() 拿到的就是异步提交的任务的返回结果
'''

# sayncio   模块实现异步IO



















