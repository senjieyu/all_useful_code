# '''多线程是否有用要看具体情况
# 单核:四个任务(O密集型\计算密集型
# 多核:四个任务(IO密集型\计算密集型
# # 计算密集型
# 每个任务都需要10s
# 单核(不用考虑了)
# 多进程:额外的消耗资源
# 多线程:介绍开销
# 多核
# 多进程:总耗时10+
# 多线程:总耗时 40+
# IO密集型
# 多核
# 多进程:相对浪费资源
# 多线程:更加节省资源'''
#
#
# # 计算密集型
# from multiprocessing import Process
# from threading import Thread
# import os,time
#
#
# def work():
#     res = 0
#     for i in range(10000000):
#         res*=i
#
# if __name__ == '__main__':
#     l=[]
#     print(os.cpu_count())
#     # 获取当前计算机CPU个数
#     start_time = time.time()
#     for i in range(32):
#         # p = Process(target=work) # 1.6276895999908447
#         t = Thread(target=work)   # 5.488417625427246
#         # p.start()
#         t.start()
#         # l.append(p)
#         l.append(t)
#     for p in l:
#         p.join()
#     print(time.time()-start_time)




# IO密集型

from multiprocessing import Process
from threading import Thread
import os,time
def work():
    time.sleep(2)


if __name__ == '__main__':
    l=[]
    print(os.cpu_count())
    # 获取当前计算机CPU个数
    start_time = time.time()
    for i in range(40):
        p = Process(target=work)   # 大 得多
        # t = Thread(target=work) #2.0807836055755615

        p.start()
#         t.start()

        l.append(p)
#         l.append(t)

    for p in l:
        p.join()
    print(time.time()-start_time)
    


'''
多进程和多线程都有各自的优势
并且我们后面在写项目的时候通常可以
多进程下面再开设多线程
这样的话既可以利用多核也可以介绍资源消耗
'''
