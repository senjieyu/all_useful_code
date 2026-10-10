import socket
import time

server = socket.socket()
server.bind(('127.0.0.1',8080))
server.listen(5)
server.setblocking(False)   # 将所有的网络阻塞变为非阻塞

r_list = []
del_list = []
while True:
    try:
        conn, addr = server.accept()
        r_list.append(conn)
    except BlockingIOError as e:
        # time.sleep(0.1)
        # print('列表的长度：',len(r_list))
        # print('做其他事')
        for conn in r_list:
            try:
                data= conn.recv(1024)  # 没有消息 报错
                if len(data) == 0: # 客户端断开链接
                    conn.close()  # 关闭conn
                    # 将无用的conn从r_list删除
                    del_list.append(conn)
                    continue
                conn.send(data.upper())
            except BlockingIOError:
                continue
            except ConnectionError:
                conn.close()
                del_list.append(conn)
        #   回收没有的链接
        for conn in del_list:
            r_list.remove(conn)
        del_list.clear()




















