#账户的输入0.退出1.登录2.注册。3.转账4.存款
u_i=input("uersername:".strip())
p_i=input('password:'.strip())
with open (r"C:\Users\senji\Desktop\python_learning\learning\ATM\log\login.txt",mode='rt',encoding="utf-8") as f:
    for x in f :    #通常在文件中使用 for 循环迭代文件对象时，默认情况下是按行取的。
        un,pw=x.strip().split(':')
        if un==u_i and pw==p_i :
            from lib.common import loggin
            def 退出():
                print("退出功能实现了")
                loggin(u_i,'退出')
                pass


            def 登录():
                from lib.common import openfile
                print("登录功能实现了")
                openfile
                loggin(u_i,'登录')
                pass

            def 注册():
                from lib.common import login
                print("注册功能实现了")
                login()
                loggin(u_i,'注册')
                pass

            def 转账():
                print("转账功能实现了")
                loggin(u_i,'转账')
                pass

            def 存款():
                print("存款功能实现了")
                loggin(u_i,'存款')
                pass

            dic = {'0': ['退出', 退出],
                    '1': ['登录', 登录],
                    '2': ['注册', 注册],
                    '3': ['转账', 转账],
                    '4': ['存款', 存款]
                    }

            def dic_use():
                for k in dic:
                    print(k, dic[k][0])

            dic_use()
            num = input('请输入服务编号：'.strip())
            break
    else: print("账号或密码错误")








