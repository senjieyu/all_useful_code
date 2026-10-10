'''
逻辑接口层
'''
from lib import common
from db import db_handle
def register_interface(username,password,balance=15000):
    user_dic=db_handle.select(username)
    if user_dic:
        return False,'用户名已存在'
    password = common.get_pwd_md5(password)
    user_dic = {
        'username': username,
        'password': password,
        'balance': 15000,
        'flow': [],
        'shop_car': {},
        'locked': False
    }
    db_handle.save(user_dic)
    return True,'注册成功'

def login_interface(username,password):
    user_dic=db_handle.select(username)
    if user_dic:
        password=common.get_pwd_md5(password)
        if user_dic['password']== password:
            return True,'登录成功'
        else:
            return False,'密码错误'
    return False,'用户不存在'

def check_bla_interface(login_user):
    res=db_handle.select(login_user)
    res1=res['balance']
    return res1

def check_flow_interface(username):
    user_dic = db_handle.select(username)
    user_flow = user_dic.get('flow')
    return user_flow


