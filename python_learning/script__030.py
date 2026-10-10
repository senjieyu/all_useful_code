# from threading import Thread,Lock
# import time
#
#
# money = 100
# mutex = Lock()
#
# def task():
#     global money
#     mutex.acquire()
#     tem = money
#     time.sleep(0.01)
#     money = tem - 1
#     mutex.release()
#
# if __name__ == '__main__':
#
#     t_list =[]
#     for i in range (100):
#         t = Thread(target=task)
#         t.start()
#         t_list.append(t)
#
#     for t in t_list:
#         t.join()
#     print(money)


from threading import Thread,Lock
import time


money = 100
mutex = Lock()

def task():
    global money
    with mutex:
        tem = money
        time.sleep(0.01)
        money = tem - 1


if __name__ == '__main__':

    t_list =[]
    for i in range (100):
        t = Thread(target=task)
        t.start()
        t_list.append(t)

    for t in t_list:
        t.join()
    print(money)
