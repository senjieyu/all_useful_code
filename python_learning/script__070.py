i_n=input ('username:')
special_characters = ['/','.', ';', ']', '[', '\\', '=']
if any(char in i_n for char in special_characters):
    print('不允许输入特殊字符')
else:
    with open( '11.txt', mode ='rt',encoding='utf-8') as f :
        for x in f:
            old_name,old_password = x.strip().split(':')
            if i_n == old_name:
                print('呆子，来晚了，已经有人用这个名字了')
                break
            else:
                i_p = input('password:')
                with open('11.txt', mode='at', encoding='utf-8') as f:
                    f.write('{}:{}\n'.format(i_n, i_p))
                    print('successful')
                    break

