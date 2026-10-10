#在调用函数的时候，直接或间接的调用到本身


#直接调用
# def f1():
#     print('sdf')
#     f1()
#
# f1()

#间接调用
# def f1():
#     f2()
# def f2():
#     f1()
# f1()


# def f1(n):
#     if n ==10:
#         return
#     print(n)
#     n=n+1
#     f1(n)
#
# f1(0)
#
#递归的本质就是循环


# #递归的两个阶段
# 回溯：一层一层调用下去
# 递推：满足某种结束条件，结束递归调用，然后一层一层返回

# def age(n):
#     if n==1:
#         return 18
#     return age(n-1)+10
# res=age(5)
# print(res)



# 应用：将列表中的值取出来
l=[1,2,[3,[4,[5,[6,[7,[8,[9,10,11]]]]]]],12]
def f(l):
    for k in l:
        if type(k) is list:
            f(k)
        else:
            print(k)

f(l)
