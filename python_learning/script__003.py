#没有实际用到，只是尝试练习


# def index(x,y):
#     print(x,y)

def deco(fun):
    def wrapper(*args,**kwargs):
        res=fun(*args,**kwargs)
        #新增功能：
        print('asss')
        #。。。。。

        return res
    return wrapper

# a=deco(index)
# a(1,5)

@deco
def index(x,y):
    print(x,y)
index(1,5)