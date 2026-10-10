# def func():
#     while True:
#         print("开始了")
#         x=yield 666
#         print("这是x的值{}".format(x))
#
# g=func()
# x=g.send(None)#第一次必须传None
# c=g.send('傻波一')
# v=g.send(1)
# g.close()#close之后将无法传值
# print(x,c,v)

def func():
    f_list=[]
    while True:
        print("开始了")
        x=yield f_list
        print("这是x的值{}".format(x))
        f_list.append(x)

g=func()
ret=g.send(None)
print(ret)

ret=g.send("时代大厦")
print(ret)