import socket
from threading import Thread
from multiprocessing import Process
'''
服务端
    1.要有固定的IP和PORT
    2.24小时不间断提供服务
    3.能够支持并发
    
养成看源码的习惯
'''
server = socket.socket()  # 括号内不加参数默认就是TCP协议
server.bind(('127.0.0.1',8080))
server.listen(5)

# 将服务的代码单独封装成一个函数
def talk (conn):
    # 通信循环
    while True:
        try:
            data = conn.recv(1024)
            # 针对mac linux 客户端
            if len(data) == 0 : break
            print(data.decode('utf-8'))
            conn.send(data.upper())
        except ConnectionError as e:
            print(e)
            break
    conn.close()

# 链接循环
while True:
    conn,addr = server.accept() # 接客
    # 叫其他人来服务客户
    t = Thread(target=talk,args=(conn,))
    t.start()






