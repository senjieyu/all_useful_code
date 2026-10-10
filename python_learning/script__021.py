# 队列Queue模块
# 队列：先进先出
# 堆栈：先进后出
# import queue
#
# # 创建一个队列
# q = queue.Queue()  # 括号内可以传数字 标示生成的队列最大可以同时存放的数据量
#
# # 往队列中存数据
# q.put(111)

'''
存取数据 存是为了更好的取
'''

from multiprocessing import Queue

q = Queue(5)

q.put(111)
q.put(2)
q.put(3)
print(q.full())  # 判断当前队列是否满了
print(q.empty())  # 判断当前队列是否kong的
q.put(4)
q.put(5)
q.put(5)  # 当队列数据放满了之后，如果还有数据要放，程序会阻塞，直到有位置让出来，。不报错

v1 = q.get(111)
v2 = q.get(111)
v3 = q.get(111)
v4 = q.get(111)
v5 = q.get(111)
v6 = q.get(111)  # q.get_nowait()  #没有数据直接报错
# q.get(timeout = 3) #没有数据之后等待三秒，后再报错
# 队列中如果已经没有数据的话 get方法会原地阻塞
print(v6, v5, v4, v3, v2, v1)
