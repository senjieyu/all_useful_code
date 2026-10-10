# 协程

# 帮助我们检测IO实现切换


import time
from gevent import spawn
from gevent import monkey
monkey.patch_all()



''' 
gevent模块本身无法检测常见的一些io操作
在使用的时候需要你额外的导入一句话
from gevent import monkey
monkey.patch_all()
又由于上面的两句话在使用gevent模块的时候是肯定要导入的
所以还支持简写
from gevent import monkey;monkey.patch_all()

'''


def heng() :
    print('哼')
    time.sleep(2)
    print('哼')
def ha() :
    print('哈')
    time.sleep(3)
    print('哈')

start_time= time.time()
g1 = spawn(heng)
g2 = spawn(ha)
g1.join()
g2.join() # 等待被检测的任务执行完毕 再往后继续执行
print(time.time() - start_time)




















