# 服务端应该满足两个特点：
# 1、一直对外提供服务
# 2、并发的服务多个客户端




from socket import  *
import subprocess
import struct
import json

server = socket(AF_INET,SOCK_STREAM)
server.bind(('127.0.0.1',8080))
server.listen(5)

#服务端应该做的两件事
#第一件事，循环的从半连接池中取出链接请求与其建立双向链接，拿到链接对象
while True:
    conn,client_addr=server.accept()

    #第二件事：拿到链接对象，与其进行通信循环
    while True:
        try:
            res=conn.recv(1024)
            if len(res) == 0:break
            obj=subprocess.Popen(res.decode('utf-8'),
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
            stderr_res=obj.stderr.read()
            stdout_res=obj.stdout.read()
            total_size = len(stderr_res)+len(stdout_res)
            # 制作头
            header_dic={
                'filename':'a.txt',
                'total_size':total_size,
                'md5':'1354464asd'
            }

            # 1.先发头信息（固定长度的bytes）：对数据描述信息
            #int  -》 固定长度的bytes     struct模块
            json_str = json.dumps(header_dic)
            json_str_bytes = json_str.encode('utf-8')

            # 2、先把头的长度发过去
            x = struct.pack('i', total_size)
            conn.send(x)


            # 3、发头信息
            conn.sendto(json_str_bytes)
            # 4.再发真实信息
            conn.send(stdout_res)
            conn.send(stderr_res)







        except Exception:
            break
    conn.close()




