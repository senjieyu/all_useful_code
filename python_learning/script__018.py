# 守护进程
# 主进程死了守护进程就死
# 在Python中，如果主进程退出，守护进程会被强制终止，而不管它们是否执行完毕。
# 这是Python解释器的一种行为，主要是为了确保程序的整体稳定性。
# 他死我就死
import time
from multiprocessing import Process


def task(name):
    print('%s总管正在活着' % name)
    time.sleep(3)
    print('%s总管正在死亡' % name)


if __name__ == '__main__':
    p = Process(target=task, args=('egon',))
    # 当形参是（name = 'asds')这种形式时可用下面的赋值形式
    # p = Process(target=task, kargs={'name': 'egon'})
    p.daemon = True  # 将进程p设置成守护进程，这句话，
    # 一定要放在p.start方法上面才有效，否则会直接报错
    p.start()
    time.sleep(2)
    print('主进程死了')