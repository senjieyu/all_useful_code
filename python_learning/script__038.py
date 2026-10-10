'''
进程:资源单位
线程:执行单位
协程:这个概念完全是程序员自己意淫出来的 根本不存在
    单线程下实现并发
    我们程序员自己再代码层面上检测我们所有的IO操作
    一旦遇到IO了 我们在代码级别完成切换
    这样给cPu的感觉是你这个程序一直在运行 没有IO
    从而提升程序的运行效率

多道技术
    切换+保存状态
    CPU两种切换
        1.程序遇到IO
        2.程序长时间占用



TCP服务端
    accept
    recv
代码如何做到
    切换+保存状态
    切换
    切换不一定是提升效率 也有可能是降低效率
    IO切   提升
    没有IO切 降低

保存状态
    保存上一次我执行的状态 下一次来接着上一次的操作继续往后执行
    yield
'''

# 验证切换是否就一定提升效率

# import time
#
#
# def func1():
#     for i in range(100000000):
#         i + 1
#
#
# def func2():
#     for i in range(100000000):
#         i + 1
#
# start_time=time.time()
# func1()
# func2()
# print(time.time() - start_time)

import time

def func1() :
    while True:
        100000000 + 1
        yield

def func2():
    g = func1()
    # 先初始化出生成器
    for i in range(100000000):
        i+1
        next(g)

start_time=time.time()
func2()
print(time.time() - start_time)

















