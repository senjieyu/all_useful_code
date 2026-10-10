# join是让主进程的代码等待子进程代码运行结束之后，再继续运行
# 不影响其他子进程的运行

#
# from multiprocessing import Process
# import time
#
#
# def task(name):
#     print('%s is running' % name)
#     time.sleep(3)
#     print('%s is over' % name)
#
#
# if __name__ == '__main__':
#     p1 = Process(target=task, args=('jason',))
#     p2 = Process(target=task, args=('egon',))
#     p3 = Process(target=task, args=('tank',))
#     p1.start()  # 告诉操作系统帮你创建一个进程 异步
#     p2.start()  # 告诉操作系统帮你创建一个进程 异步
#     p3.start()  # 告诉操作系统帮你创建一个进程 异步
#     # p.join()  # 主进程等待子进程p运行结束之后再继续往后执行
#     print('zhu')


from multiprocessing import Process
import time


def task(name, n):
    print('%s is running' % name)
    time.sleep(n)
    print('%s is over' % name)


if __name__ == '__main__':
    # 1
    p1 = Process(target=task, args=('jason', 1))
    p2 = Process(target=task, args=('egon', 2))
    p3 = Process(target=task, args=('tank', 3))
    start_time = time.time()
    p1.start()  # 告诉操作系统帮你创建一个进程 异步
    p2.start()  # 告诉操作系统帮你创建一个进程 异步
    p3.start()  # 告诉操作系统帮你创建一个进程 异步
    p1.join()  # 主进程等待子进程p1运行结束之后再继续往后执行
    p2.join()  # 主进程等待子进程p2运行结束之后再继续往后执行
    p3.join()  # 主进程等待子进程p3运行结束之后再继续往后执行
    print('zhu', time.time() - start_time)

    #   2
    # start_time = time.time()
    # for i in range(1,4):
    #     p = Process(target=task,args=('子进程%s' %i, i))
    #     p.start()
    #     p.join()
    # print('zhu', time.time() - start_time)

    #  3
    # start_time = time.time()
    # p_list = []
    # for i in range(1,4):
    #     p = Process(target=task,args=('子进程%s' %i, i))
    #     p.start()
    #     p_list.append(p)
    # # 确保先启动三个进程
    # for p in p_list:
    #     p.join()
    # print('zhu', time.time() - start_time)

# 进程之间数据相互隔离
