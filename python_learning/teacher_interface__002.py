from lib import common
from db import db_handler
def register_interface(username,password):
    user_dic=db_handler.select(username,'老师')
    if user_dic:
        return False,'用户名已存在'
    password = common.get_pwd_md5(password)
    user_dic = {
        'username': username,
        'password': password,
        '班级':[],
        '上课班级':'待定',
        '学生列表':[]
    }
    db_handler.save(user_dic,'老师')
    return True,'注册成功'


def login_interface(username,password):
    user_dic=db_handler.select(username,'老师')
    if user_dic:
        password=common.get_pwd_md5(password)
        if user_dic['password']== password:
            return True,'登录成功'
        else:
            return False,'密码错误'
    return False,'用户不存在'
