
#记录日志功能
def loggin(user,a):
    import time
    with open (r'C:\Users\senji\Desktop\python_learning\learning\ATM\log\user.log',mode='at',encoding='utf-8') as f:
        f.write('用户{}于{}进行{}操作了\n'.format(user,time.strftime('%Y年%m月%d日%H时%M分%S秒'),a))

#登录操作
def openfile():
    u_i=input("uersername:")
    p_i=input('password:')
    with open (r"C:\Users\senji\Desktop\python_learning\learning\ATM\log\login.txt",mode='rt',encoding="utf-8") as f:
        for x in f :    #通常在文件中使用 for 循环迭代文件对象时，默认情况下是按行取的。
            un,pw=x.strip().split(':')
            if un==u_i and pw==p_i :
                print('succeedful')
                break
        else: print("账号或密码错误")

#注册操作
def login():
    i_n=input ('username:')
    special_characters = ['/', '.', ';', ']', '[', '\\', '=']
    if any(char in i_n for char in special_characters):
        print('不允许输入特殊字符')
    else:
        with open( r'C:\Users\senji\Desktop\python_learning\learning\ATM\log\login.txt',mode ='rt',encoding='utf-8') as f :
            for x in f:
                old_name,old_password = x.strip().split(':')
                if i_n == old_name:
                    print('呆子，来晚了，已经有人用这个名字了')
                    break
                else:
                    i_p = input('password:')
                    with open(r'C:\Users\senji\Desktop\python_learning\learning\ATM\log\login.txt', mode='at', encoding='utf-8') as f:
                        f.write('{}:{}\n'.format(i_n, i_p))
                        print('successful')
                        break