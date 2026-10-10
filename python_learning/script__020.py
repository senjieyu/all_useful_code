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
# """
# Windows操作系统下，创建进程一定要在main内创建
# 因为Windows下创建进程类似于模块导入的方式
# 会从上往下依次执行代码
# """
# if __name__ == '__main__':
#     # 创建进程的代码
#     # 1.创建一个对象
#     p = Process(target=task, args=('jason',))
#     # 容器类型哪怕里面只有一个元素 一定要用逗号隔开
#     # 2.开启进程
#     p.start()  # 告诉操作系统帮你创建一个进程
#     print('zhu')
#
# # 总结 进程就是在内存中申请一块内存空间将需要运行的代码丢进去
# # 一个进程对应在内存中就是一块独立的内存空间
# # 多个进程对应在内存中就是多块独立的内存空间
# # 进程与进程之间数据默认情况下是无法直接交互，如果想交互可以借助第三方工具、模块
#
import time
# 第二种方式 类的继承式

from multiprocessing import Process


class MyProcess(Process):
    def run(self):   #必须叫run
        print('11111110')
        time.sleep(1)
        print('222222222')

if __name__ == '__main__':
    p = MyProcess()
    p.start()
    print('zhu')








