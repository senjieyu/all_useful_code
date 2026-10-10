# class A(object):
#     def test(self):
#         print('from A')
#
#
# class B(A):
#     def test(self):
#         print('from B')
#
#
# class C(A):
#     def test(self):
#         print('from C')
#
#
# class D(B,C):
#     pass
#
#
# obj = D()
# obj.test() # 结果为：from B
# print(D.mro())  # 新式类内置了mro方法可以查看线性列表的内容，经典类没有该内置该方法
#
# # python会在MRO列表上从左到右开始查找基类,直到找到第一个匹配这个属性的类为止。 而这个MRO列表的构造是通过一个C3线性化算法来实现的。我们不去深究这个算法的数学原理,它实际上就是合并所有父类的MRO列表并遵循如下三条准则:
# #
# # 1.子类会先于父类被检查
# # 2.多个父类会根据它们在列表中的顺序被检查
# # 3.如果对下一个类存在两个合法的选择,选择第一个父类
# # 所以obj.test()的查找顺序是，先从对象obj本身的属性里找方法test，没有找到，则参照属性查找的发起者(即obj)所处类D的MRO列表来依次检索，首先在类D中未找到，然后再B中找到方法test




# 参照下述代码，多继承结构为非菱形结构，此时，
# 会按照先找B这一条分支，然后再找C这一条分支，
# 最后找D这一条分支的顺序直到找到我们想要的属性
class E:
    def test(self):
        print('from E')


class F:
    def test(self):
        print('from F')


class B(E):
    def test(self):
        print('from B')


class C(F):
    def test(self):
        print('from C')


class D:
    def test(self):
        print('from D')


class A(B, C, D):
    # def test(self):
    #     print('from A')
    pass


print(A.mro())
'''
[<class '__main__.A'>, <class '__main__.B'>, <class '__main__.E'>, <class '__main__.C'>, <class '__main__.F'>, <class '__main__.D'>, <class 'object'>]
'''

obj = A()
obj.test() # 结果为：from B
# 可依次注释上述类中的方法test来进行验证
