'''
这里研究的IO模型都是针对网络IO的
我们这里研究的Io模型都是针对网络Io的stevens在文章中一共比较了五种IO Model:
* blocking Io
阻塞Io
*nonblocking Io
非阻塞Io
*1o multiplexing
Io多路复用
*signal driven Io
信号驱动Io
* asynchronous Io异步IO
由signal driven Io (信号驱动Io)在实际中并不常用，所以主要介绍其余四种IOModel。
'''


# 1)等待数据准备 (waiting for the data to be ready)
# 2) 将数据内核拷贝到进程中(Copying the data from the kernel to the process

# 同步异步
# 阻塞非阻塞
# 常见的网络阻塞状态:
#     accept
#     recv
#     recvfrom
#     send虽然它也有io行为 但是不在我们的考虑范围

# 阻塞 IO   之前写的都是
# 在服务端开设多进程或者多线程 进程池线程池 其实还是没有解决IO问题
# 该等的地方还是得等 没有规避
# 只不过多个人等待的彼此互不干扰


# 非阻塞 IO

'''
虽然非阻塞IO给你的感觉非常的牛逼
但是该模型会
长时间占用着CPu并且不干活 让cPu不停的空转
我们实际应用中也不会考虑使用非阻塞IO模型
'''

# IO多路复用
'''
当监管的对象只有一个的时候 其实IO多路复用连阻塞IO都比比不上!!!
但是IO多路复用可以一次性监管很多个对象
监管机制是操作系统本身就有的 如果你想要用该监管机制(select)
需要你导入对应的select模块
'''


# 异步 IO

# 见 https://zhuanlan.zhihu.com/p/112185254









