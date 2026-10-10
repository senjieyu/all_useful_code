# 生产者：生产、制造东西的
# 消费者：消费、处理东西的
# 该模型除了上述两个之外还需要一个媒介
# 生产者+消息队列+消费者
from multiprocessing import Process, Queue
import time
import random


def producer(name, food, q):
    for i in range(4):
        data = '%s生产了%s%s' % (name, food, i)
        time.sleep(random.randint(1, 3))
        print(data)
        # 将数据放入队列中
        q.put(data)


def consumer(name, q):
    # 消费者胃口很大，光盘行动
    while True:
        food = q.get()  # 没有数据就会卡住
        # 判断当前是否有结束的标示
        if food is None:break

        time.sleep(random.randint(1, 3))
        print('%s吃了%s' % (name, food))


if __name__ == '__main__':
    q = Queue()
    p1 = Process(target=producer, args=('大厨', '包子', q))
    p2 = Process(target=producer, args=('jack', 'soup', q))
    c1 = Process(target=consumer, args=('ming', q))
    c2 = Process(target=consumer, args=('xxxx', q))
    p1.start()
    p2.start()
    c1.start()
    c2.start()

    # 等待生产者生产完毕
    p1.join()
    p2.join()
    # 等待生产者生产完毕之后往队列中添加特定的结束符号
    q.put(None)  # 肯定在所有生产者生产的数据的末尾
    q.put(None)  # 肯定在所有生产者生产的数据的末尾
