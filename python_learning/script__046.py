
#
# lambda x,y:x+y
# print(lambda x,y:x+y)
#
#
# # 调用方式1.
# # 内存地址加括号
# res=(lambda x,y:x+y)(1,2)
# print(res)

# # 方式二
# func=lambda x,y:x+y
# print(func(1,2))


#匿名函数，没有名字，不赋值的话，引用计数为0，运行下一行代码时便会销毁
#节省内存空间，用于临时调用一次的场景：
#更多的是将匿名函数与其他函数配合使用

# l={'asd':465,
#    'asf':4866,
#    'wad':1599,
#    'adsfe':499}
# #需求：找出薪资最高的
# res=max(l,key=lambda k:l[k])
# #max函数默认是按照字典的key的排序
# #res=max(l)#默认的排序得不到正确的结果
# print(res)
# res=min(l,key=lambda k:l[k])
# print(res)
# res=sorted(l,key=lambda k:l[k])#规定
# print(res)


# l=['asd','erh','sdab','asega']
# l1=[name+'+dsb' for name in l]
# print(l1 )

# l=['asd','erh','sdab','asega']
# l1=[name+'+dsb' for name in l]
# res=map(lambda name:name+'+dsb',l1)#了解
# print(res)#为生成器

#上面同此
# l=['asd','erh','sdab','asega']
# l1=(name+'+dsb' for name in l)#生成器生成式
# print(l1 )


