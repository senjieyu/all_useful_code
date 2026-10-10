from multiprocessing import Process, Lock
import json
import time
import random

'''
多个进程操作同一份数据时，会出现数据错乱的问题
针对上述问题，解决方式就是加锁处理：
将并行变成串行，牺牲效率但是保证了数据安全
'''


# 查票
def search(i):
    # 文件操作读取票数
    with open('data', 'r', encoding='utf-8') as f:
        dic = json.load(f)
    print('用户%s查询余票%s' % (i, dic.get('ticket_num')))
    # 字典取值不要用[]的形式 推荐使用get  你写的代码打死都不能报错


# 买票 1.先查  2.后买
def buy(i):
    # 先查票
    with open('data', 'r', encoding='utf-8') as f:
        dic = json.load(f)
    # 模拟网络延迟
    time.sleep(random.randint(1, 3))
    # 判断当前是否有票
    if dic.get('ticket_num') > 0:
        # 修改数据库 买票
        dic['ticket_num'] -= 1
        # 写入数据库
        with open('data', 'w', encoding='utf-8') as f:
            json.dump(dic, f)

        print('用户%s买票成功' % i)
    else:
        print('用户%s买票失败' % i)


def run(i,mutex):
    search(i)
    # 给买票环节枷锁处理
    # 抢锁
    mutex.acquire()
    buy(i)
    mutex.release()

if __name__ == '__main__':
    # 再主进程中生成一把锁，让所有的子进程抢，谁先抢到谁先买票
    mutex = Lock()
    for i in range(1, 11):
        p = Process(target=run, args=(i, mutex))
        p.start()


# 1.锁不要轻易使用，容易造成死锁现象（我们写代码一般不会用，
# 都是内部封装好的）
# 2.锁只在处理数据的部分加来保证数据安全（只在争取数据的环节加锁处理即可）


