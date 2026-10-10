
i_io=input('请输入你的文件的绝对路径：')
i_out=input('请输入你的文件的绝对路径：')
with open (r'{}'.format(i_io),mode='rt',encoding='utf-8') as f1,\
        open (r'{}'.format(i_out),mode='wt',encoding='utf-8') as f2:
    read_old = f1.read()
    new_words = f2.write(read_old)

