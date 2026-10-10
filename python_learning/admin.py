from interface import  admin_interface
from lib import common
from interface import common_interface


admin_info={
    'user':None
}
def register():
    while True:
        username= input('请输入用户名：').strip()
        password=input('请输入密码：').strip()
        re_password = input('请确认密码：').strip()
        if password== re_password:
            flag,msg=admin_interface.admin_register_interface(
                username,password
            )
            if flag:
                print(msg)
                break
            else:
                print(msg)

        else:
            print('两次密码不一致，请重新数入')



def login():
    while True:
        username = input('请输入用户名：').strip()
        password = input('请输入密码：').strip()

        flag, msg = common_interface.login_interface(
            username, password,user_type='admin'
        )
        if flag:
            print(msg)
            admin_info['user'] = username
            break
        else:
            print(msg)




@common.auth('admin')
def creat_school():
    while True:
        school_name = input('请输入学校名称：').strip()
        school_addr = input('请输入学校地址：').strip()
        flag,msg = admin_interface.create_school_interface(
            school_name,school_addr,admin_info['user']
        )
        if flag:
            print(msg)
            break
        else:
            print(msg)





@common.auth('admin')
def creat_course():
    # 让管理员先选择学校，
    # 选择学校后，再输入课程名称
    # 调用创建课程接口，让管理员去创建课程
    while True:
        flag,school_list_or_msr = common_interface.get_all_school_interface()
        if not flag:
            print(school_list_or_msr)
            break
        for index,school_name in enumerate(school_list_or_msr):
            print(index,school_name)
        choice = input('请输入学校编号：').strip()
        if not choice.isdigit():
            print('请输入数字')
            continue
        choice = int(choice)
        if choice not in range(len(school_list_or_msr)):
            print('请输入正确编号！')
            continue

        school_name = school_list_or_msr[choice]
        course_name = input('请输入需要创建的可成名称：').strip()
        flag,msg = admin_interface.create_course_interface(
            school_name,course_name,admin_info['user']
        )
        if flag:
            print(msg)
            break
        else:
            print(msg)


@common.auth('admin')
def creat_teacher():
    while True:
        teacher_name = input('请输入老师的名字：').strip()
        flag,msg = admin_interface.create_teacher_interface(
            teacher_name,admin_info['user']
        )
        if flag:
            print(msg)
            break
        else:
            print(msg)





func_dict={
    '1':register,
    '2':login,
    '3':creat_school,
    '4':creat_course,
    '5':creat_teacher

}






def admin_view():
    while True:
        print('''
        管理员
        1.注册
        2.登录
        3.创建学校
        4.创建课程
        5.创建讲师
        ''')

        choice = input('请输入功能编号：').strip()
        if choice == 'q':
            break

        if choice not in func_dict:
            print('输入有误，请重新输入')
            continue

        func_dict.get(choice)()
















