from socket import *
client = socket(AF_INET,SOCK_STREAM)
client.connect(('127.0.0.1',8080))
import struct
import json

while True:
    msg=input('请输入命令：').strip()
    if len(msg) == 0 :continue
    client.send(msg.encode('utf-8'))
    # 解决粘包问题思路：
    #一、先收固定长度的头：解析出数据的描述信息，包括数据的总大小
    header=client.recv(4)

    total_size=struct.unpack('i',header)[0]

    # 二、根据解析出的描述信息，接收真实的数据
    #2、recv_size=0，循环接收，每接受一次，recv_size+=接收的长度
    #3、直到recv_size=total_size,结束循环

    # 接收端
    # 1、先收4个字节，从中提取接下来要收的头的长度
    x = client.recv(4)
    header_len = struct.unpack('i',x)[0]


    # 2.接收头，并解析
    json_str_bytes = client.recv(header_len)
    json_str = json_str_bytes.decode('utf-8')
    header_dic = json.loads(json_str)
    print(header_dic)
    total_sice= header_dic['total_size']

    #3.接收真实的数据
    recv_size =0
    while recv_size<total_size:
        recv_data=client.recv(1024)#本次接收，最大接收1024个字节
        recv_size += len(recv_data)
        print(recv_data.decode('gbk'),end='')  # 强调：此处用的系统编码格式
    else:
        print()




#粘包问题出现的原因
#1.tcp是流式协议，数据像水流一样粘在一起，没有任何边界区分
#2.收数据没收干净，有残留，就会下一次结果混淆在一起
#解决的核心法门就是：每次都收干净，不要任何残留



