

#装饰器模板总结
def func():
    def wrapper(*args,**kwargs):
        pass
        ret=func(*args,**kwargs)
        pass
        return ret
    return wrapper


