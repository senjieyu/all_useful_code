
login_teacher=False
login_student=False
login_ad=False
from interface import student_interface
from interface import teacher_interface
from interface import Administrator_interface
# 管理员操作

from lib.common import ad_login_auth

@ad_login_auth
def 系统初始化():
    chose=input('请确认要进行初始化操作，确认请按 1：')
    if chose=='1' :
        chose = input('请再次确认要进行初始化操作，确认请按 1：')
        if chose == '1':
            chose = input('请最终确认要进行初始化操作，确认请按 2：')
            if chose == '2':
                Administrator_interface.renew()
                print('操作成功')
            else:
                print('操作取消')
        else:
            print('操作取消')
    else:
        print('操作取消')

def 管理员登录():
    while True:
        username = input('请输入用户名:').strip()
        password = input('请输入密码：').strip()
        flag, res = Administrator_interface.login_interface(username, password)
        if flag:
            print(res)
            global login_ad
            login_ad = username
            break
        else:
            print(res)
            break

@ad_login_auth
def 创建班级():
    while True:
        chose = input('创建班级请按1，删除班级请按2，查看班级请按3，退出请按任意键')
        if chose == '1':
            ad_creat_class = Administrator_interface.CreatedClass()
            ad_creat_class.show_class_and_school()
            school_name=input("请输入您要创建的校区：").strip()
            flag1,res1=ad_creat_class.chose_school(school_name)
            if flag1:
                print(res1)
                flag = True
                while flag:
                    class_name = input('请输入您要创建的班级名称：').strip()
                    flag, res = ad_creat_class.creat_class(school_name,class_name)
                    if flag:
                        print(res)
                        ad_creat_class.show_class_and_school()
                        flag2 = input('不再添加请按任意键，继续请按1：')
                        if flag2 =='1':
                            continue
                        else:
                            break
                    else:
                        print(res)
                        break
                    break
            else:
                print(res1)
            continue

        elif chose == '2':
            ad_creat_class = Administrator_interface.CreatedClass()
            ad_creat_class.show_class_and_school()
            school_name = input("请输入您要创建的校区：").strip()
            class_name = input('请输入您要删除的班级名称：').strip()
            flag,res = ad_creat_class.del_class(school_name,class_name)
            if flag:
                print(res)
            else:
                print(res)
            ad_creat_class.show_class_and_school()
            continue
        elif chose == '3':
            ad_creat_class.show_class_and_school()
            continue
        else:
            break




@ad_login_auth
def 管理学校():
    while True:
        ad_set = Administrator_interface.AdSet()
        chose = input('创建学校请按1，删除学校请按2，查看学校请按3，退出请按任意键')
        if chose == '1':
            school_name = input('请输入您要创建的学校名称：').strip()
            ad_set.add_school(school_name)
            ad_set.show_school()
            continue
        elif chose == '2':
            school_name = input('请输入您要删除的学校名称：').strip()
            ad_set.remove_school(school_name)
            ad_set.show_school()
            continue
        elif chose == '3':
            ad_set.show_school()
        else:
            break


@ad_login_auth
def 创建讲师():



    pass

@ad_login_auth
def 创建课程():
    print('ok')

ad_dict = {
    '0': 管理员登录,
    '1': 管理学校,
    '2': 创建班级,
    '3': 创建讲师,
    '4': 创建课程,
    '5':5,
    '666':系统初始化
}


def ad_view():
    while True:
        print("""
                    0.管理员登录
                    1.管理学校
                    2.管理班级
                    3.创建讲师
                    4.创建课程
                    5.返回上一界面
                    666.系统初始化
        """)
        ad_number = input('请选择你要进行的业务：')
        if ad_number in ad_dict:
            if ad_number =='5':
                break
            else:
                ad_dict[ad_number]()
        else:
            print("请输入正确的数字编号")
#======================================================================================================
#======================================================================================================
#======================================================================================================
#======================================================================================================
#======================================================================================================
#======================================================================================================

#老师操作：
from lib.common import tea_login_auth
def 老师注册():
    while True:
        username = input('请输入用户名:').strip()
        password = input('请输入密码：').strip()
        re_password = input('请确认密码：').strip()
        if password == re_password:
            flag, msg = teacher_interface.register_interface(username, password)
        else:
            print('请确保密码前后一致')
            continue
        if flag:
            print(msg)
            break
        else:
            print(msg)
            break

def 老师登录() :
    while True:
        username = input('请输入用户名:').strip()
        password = input('请输入密码：').strip()
        flag, res = teacher_interface.login_interface(username, password)
        if flag:
            print(res)
            global login_teacher
            login_teacher = username
            break
        else:
            print(res)
            break

@tea_login_auth
def 选择班级():
    print('ok')

@tea_login_auth
def 查看学员列表():
    print('ok')



tea_dict={
                    '0':老师注册,
                    '1':老师登录,
                    '2':选择班级,
                    '3':查看学员列表,
                    '4':4
}

def teacher_view():
    while True:
        print("""
                    0.注册
                    1.登录
                    2.选择班级
                    3.查看学员列表
                    4.返回上一界面
        """)
        tea_number=input('请选择你要进行的业务：')
        if tea_number in stu_dict:
            if tea_number =='4':
                break
            else:
                tea_dict[tea_number]()
        else:
            print("请输入正确的数字编号")


#======================================================================================================
#======================================================================================================
#======================================================================================================
#======================================================================================================
#======================================================================================================
# =====================================================================================================
from lib.common import stu_login_auth
#学生操作：
def 学生注册() :
    while True:
        username = input('请输入用户名:').strip()
        password = input('请输入密码：').strip()
        re_password = input('请确认密码：').strip()
        if password == re_password:
            flag, msg = student_interface.register_interface(username, password)
        else:
            print('请确保密码前后一致')
            continue
        if flag:
            print(msg)
            break
        else:
            print(msg)
            break


def 学生登录():
    while True:
        username = input('请输入用户名:').strip()
        password = input('请输入密码：').strip()
        flag, res = student_interface.login_interface(username, password)
        if flag:
            print(res)
            global login_student
            login_student = username
            break
        else:
            print(res)
            break

@stu_login_auth
def 学生交学费():
    print('ok')

@stu_login_auth
def 学生选择班级():
    print('ok')

stu_dict={
                    '0':学生注册,
                    '1':学生登录,
                    '2':学生交学费,
                    '3':学生选择班级,
                    '4':4
}

def student_view():
    while True:
        print("""
                    0.注册
                    1.登录
                    2.交学费
                    3.选择班级
                    4.返回上一界面
        """)
        stu_number=input('请选择你要进行的业务：')
        if stu_number in stu_dict:
            if stu_number =='4':
                break
            else:
                stu_dict[stu_number]()
        else:
            print("请输入正确的数字编号")


# =====================================================================================================
#======================================================================================================
#======================================================================================================
#======================================================================================================
#======================================================================================================



def run():
    while True:
        print("""
                    1.老师
                    2.学生
                    3.管理员
        """)
        number=input("请选择您的身份:")
        if number in ['1','2','3']:
            if number == '1':
                teacher_view()
            elif number == '2':
                student_view()
            elif number == '3':
                ad_view()
        else:
            print("请按要求属于对应数字")