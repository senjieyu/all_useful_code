#应用原来的功能之前先输密码，账号、密码在库文件中再进行

import time


def password(func):
    def wrapper(*args,**kwargs):
        pass_words=input('请输入账号：').strip()
        code=input("请输入密码：").strip()
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

    return wrapper


@password
def index(x,y):
    print('输出的是：{}'.format(x))
    time.sleep(0.1)
    print('输出的是：{}'.format(y))
    return 123

f1=index(12,56)
print(f1)