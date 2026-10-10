'''
用户视图层
'''
import os.path
from interface import user_interface
#设置全局变量，记录用户是否登录
login_user= None
from lib import common
from interface import bank_interface
# 1.注册功能
'''面条版'''
'''
def register():
    while True:
        username=input('请输入用户名:').strip()
        password=input('请输入密码：').strip()
        re_password=input('请确认密码：').strip()

        if password==re_password :
            user_dic={
                'username':username,
                'password':password,
                'balance':15000,
                'flow':[],
                'shop_car':{},
                'locked': False
            }
            import json
            from conf.settings import BASE_PATH
            use_path=os.path.join(BASE_PATH,f'{username}.json')

            # 查看用户是否存在
            if os.path.exists(use_path):
                print("用户名已存在请重新输入！")
                with open(use_path,'r',encoding='utf-8') as f:
                user_dic=json.load(f)
                continue
            #写入用户
            with open(use_path,'w',encoding='utf-8') as f:
                json.dump(user_dic,f,ensure_ascii=False)#ensure_ascii =Flase 存储的中文可以显示
                break
'''
'''拆分版'''
def register():
    while True:
        username=input('请输入用户名:').strip()
        password=input('请输入密码：').strip()
        re_password=input('请确认密码：').strip()
        if password==re_password :
            flag,msg=user_interface.register_interface(username,password)
        if flag:
            print(msg)
            break
        else:
            print(msg)
            break

# 2.登录功能
def login():
    while True:
        username = input('请输入用户名:').strip()
        password = input('请输入密码：').strip()
        flag,res=user_interface.login_interface(username,password)
        if flag:
            print(res)
            global login_user
            login_user=username
            break
        else:
            print(res)
            break

# 3.查看余额
@common.login_auth
def check_blance() :
    balance=user_interface.check_bla_interface(login_user)
    print('账户余额为：{}'.format(balance))

# 4.提现功能
@common.login_auth
def withdraw():
    while True:
        input_money=input('请输入提现金额：').strip()
        if not input_money.isdigit():
            print('只能输入数字，请重新输入')
            continue
        input_money=int(input_money)
        flag,msg=bank_interface.withdraw_interface(login_user,input_money)
        if flag:
            print(msg)
            break
        else:
            print(msg)

# 5.还款功能
@common.login_auth
def repay():
    while True:
        money=input('请输入存款金额：').strip()
        flag,msg=bank_interface.bank_repay(login_user,money)
        if flag:
            print(msg)
            break
        else:
            print(msg)

# 6.转账功能
@common.login_auth
def transfer():
    while True:
        transfer_user=input('请输入转入账户名:')
        re_sure=input('请确认转入账户：')
        if transfer_user==re_sure:
            while True:
                money = input('请输入转入金额：')
                if money.isdigit() and money > '0':
                    flag,msg=bank_interface.bank_interface_transfer(login_user,transfer_user,money)
                    if flag:
                        print(msg)
                        break
                    else:
                        print(msg)
                        break
                else:
                    print('请输入正确的金额')
            break
        else:
            print("账户名输入不一致")


# 7.查看流水
@common.login_auth
def check_flow():
    msg=user_interface.check_flow_interface(login_user)
    print(msg)


# 8.购物功能
@common.login_auth
def shopping():
    pass


# 9.查看购物车
@common.login_auth
def check_shop_car():
    pass


# 10.管理员功能
def admin():
    pass

#创建函数功能字典
func_dic={
    '1':register,
    '2':login,
    '3':check_blance,
    '4':withdraw,
    '5':repay,
    '6':transfer,
    '7':check_flow,
    '8':shopping,
    '9':check_shop_car,
    '10':admin
}

def run():
    while True:
        print('''                 
        
        =======ATM+购物车=========
                 1.注册功能
                 2.登录功能
                 3.查看余额
                 4.提现功能
                 5.还款功能
                 6.转账功能
                 7.查看流水
                 8.购物功能
                 9.查看购物车
                 10.管理员功能
        ===========end===========
                 ''')
        choice=input('请输功能编号：').strip()

        if choice not in func_dic:
            print('请输入正确的功能编号！')
            continue
        func_dic[choice]()  #gunc_dic.get(choice)()


