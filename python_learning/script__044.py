


def 有参装饰器(x):
    def decorator(func):
        def wrapper(*args, **kwargs):
            # 在函数调用前执行一些操作，使用参数 x
            print(f"装饰器参数 x 的值为: {x}")
            ret = func(*args, **kwargs)
            # 在函数调用后执行一些操作
            return ret
        return wrapper
    return decorator

@有参装饰器(5)  # 使用装饰器，并传递参数 5
def 被装饰对象():
    print("被装饰的函数")

# 调用被装饰的函数
被装饰对象()
