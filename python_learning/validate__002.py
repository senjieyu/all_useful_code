import time


def password(func):
    def wrapper(*args,**kwargs):
        pass_words=input('请输入账号：').strip()
        code=input("请输入密码：").strip()
        bd_type=input('请输入账号库信息：')
        #1.验证来自于文件
        if bd_type=='file':
            with open('11.txt',mode='rt',encoding='utf-8') as f1:
                line = f1.readline()
                while line:
                    pp, cc = line.strip().split(':')
                    # 在这里处理 pp 和 cc
                    line = f1.readline()
                    if pp==pass_words and cc==code :
                        ret=func(*args,**kwargs)
                        return ret
                        break
                    elif line:
                        continue
                    else :
                        print('账号或密码错误')
        #2. 验证来自于数据库
        elif bd_type=='mysql':
            pass
            ret = func(*args, **kwargs)
            return ret

        #
        elif bd_type=='Ldap':
            pass
            ret = func(*args, **kwargs)
            return ret


        else:
            pass





    return wrapper


@password
def index(x,y):
    print('输出的是：{}'.format(x))
    time.sleep(0.1)
    print('输出的是：{}'.format(y))
    return 123

f1=index(12,56)
print(f1)