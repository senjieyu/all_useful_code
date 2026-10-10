
#密码加盐
import hashlib
def get_pwd_md5(password):
    md5_obj=hashlib.md5()
    md5_obj.update(password.encode('utf-8'))
    salt="asafAsffag"
    md5_obj.update(salt.encode('utf-8'))
    return md5_obj.hexdigest()


# 登录功能的装饰器
def login_auth(func):
    from core import src
    def inner(*args,**kwargs):
        if src.login_user:
            res = func(*args, **kwargs)
            return res
        else:
            print("请先登录")
            src.login()
    return inner