from multiprocessing import Process, current_process
import time
import os

'''一台计算机上面运行着很多进程，那么计算机是如何区分管理这些进程服务端的呢？
计算机会给每一个运行的进程分配一个PID号
如何查看
cmd 输入tasklist即可查看
tasklist|findstr PID号   查看具体的进程


'''


##方法1
# def task():
#     # 查看当前进程的进程号  current_process().pid
#     print('%s is running' %current_process().pid)
#     time.sleep(3)
#
# if __name__ == '__main__':
#     p = Process(target=task)
#     p.start()
#     print('zhu')


##方法2
# def task():
#     # 查看当前进程的进程号  current_process().pid
#     print('%s is running' % os.getpid())
#     print('子进程的主进程号%s' % os.getppid())
#     time.sleep(3)
#
#
# if __name__ == '__main__':
#     p = Process(target=task)
#     p.start()
#     print('zhu', os.getpid())
#     print('zhuzhu', os.getppid())#获取父进程的pid


def task():
    # 查看当前进程的进程号  current_process().pid
    print('%s is running' % os.getpid())
    # print('子进程的主进程号%s' % os.getppid())
    time.sleep(3)


if __name__ == '__main__':
    p = Process(target=task)
    p.start()
    p.terminate()  # 杀死当前进程
    #  是告诉操作系统帮你杀死当前进程，但是需要一定的时间
    # 而代码的运行速度很快
    time.sleep(1)
    print(p.is_alive())  # 判断当前进程是否存活
    print('zhu')
    # print('zhu', os.getpid())
    # print('zhuzhu', os.getppid())#获取父进程的pid


