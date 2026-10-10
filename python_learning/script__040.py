# 但凡对象的内置函数中含有__iter__的都称为可迭代对象
# 调用可迭代对象的__iter__方法会将其转换成迭代器对象



# l=[1,2,3,4,5,6,7,8,9]
# p=l.__iter__()
# print(p.__next__())
# print(p.__next__())
# print(p.__next__())
# print(p.__next__())
# print(p.__next__())
# print(p.__next__())
# print(p.__next__())
# print(p.__next__())

q={1,2,3,4,5,6,7,8,9}
s=q.__iter__()
while True:
    try:
        print(s.__next__())
    except StopIteration:
        break





# 可迭代对象与迭代器对象详解
# 可迭代对象（“可以转换成迭代器的对象”）：内置有__iter__方法对象
# 迭代器对象：内置有__next__方法且内置有__iter__方法的对象
#    迭代器对象.__next__():得到迭代器的下一个值
#    迭代器对象.__iter__():得到迭代器本身，说白了调了跟没调一样

#优点：无论是有索引还是没有索引都可用   2.省内存
 

#for 循环的工作原理