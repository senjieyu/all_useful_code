def open():
    u_i=input("uersername:")
    p_i=input('password:')
    with open (r"11.txt",mode='rt',encoding="utf-8") as f:
        for x in f :    #通常在文件中使用 for 循环迭代文件对象时，默认情况下是按行取的。
            un,pw=x.strip().split(':')
            if un==u_i and pw==p_i :
                print('succeedful')
                break
        else: print("账号或密码错误")



