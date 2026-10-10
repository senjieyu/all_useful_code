import hashlib
def get_pwd_md5(password):
    md5_obj=hashlib.md5()
    md5_obj.update(password.encode('utf-8'))
    salt="asafAsffag"
    md5_obj.update(salt.encode('utf-8'))
    return md5_obj.hexdigest()


def tea_login_auth(func):
    from core import src
    def inner(*args,**kwargs):
        if src.login_teacher:
            res = func(*args, **kwargs)
            return res
        else:
            print("请先登录")
            src.老师登录()
    return inner

def stu_login_auth(func):
    from core import src
    def inner(*args,**kwargs):
        if src.login_student:
            res = func(*args, **kwargs)
            return res
        else:
            print("请先登录")
            src.学生登录()
    return inner

def ad_login_auth(func):
    from core import src
    def inner(*args,**kwargs):
        if src.login_ad:
            res = func(*args, **kwargs)
            return res
        else:
            print("请先登录")
            src.管理员登录()
    return inner