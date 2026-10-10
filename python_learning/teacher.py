from lib import common
from interface import common_interface
from interface import teacher_interface
teacher_info={
    'user':None
}

def login():
    while True:
        username = input('请输入用户名：').strip()
        password = input('请输入密码：').strip()

        flag, msg = common_interface.login_interface(
            username, password, user_type='student'
        )
        if flag:
            print(msg)
            teacher_info['user'] = username
            break
        else:
            print(msg)

@common.auth('teacher')
def check_course():
    flag, course_list = teacher_interface.check_course_interface(
        teacher_info.get('user')
    )
    if flag:
        print(course_list)
    else:
        print(course_list)


@common.auth('teacher')
def choose_course():
    # 1.先打印所有学校，并选择
    # 2.从选择的学校中获取所有课程
    while True:
        flag, school_list = common_interface.get_all_school_interface()
        if not flag:
            print(school_list)
            break

        for index, school_name in enumerate(school_list):
            print(f'编号：{index},学校名：{school_name}')

        choice = input('请输入选择的学校编号：').strip()
        if not choice.isdigit():
            print('输入有误')
            continue
        choice = int(choice)

        if choice not in range(len(school_list)):
            print('输入编号有误！')
            continue

        school_name = school_list[choice]

        flag1,course_list = common_interface.get_course_in_school_interface(
            school_name
        )
        if not flag1:
            print(course_list)
            break

        for index2, course_name in enumerate(course_list):
            print(f'编号：{index2},学校名：{course_name}')

        choice2 = input('请输入选择的课程编号：').strip()
        if not choice2.isdigit():
            print('输入有误')
            continue
        choice2 = int(choice2)

        if choice2 not in range(len(course_list)):
            print('输入编号有误！')
            continue

        course_name = course_list[choice2]
        flag2,msg = teacher_interface.add_course_interface(
            course_name,teacher_info['user']
        )
        if flag2:
            print(msg)
            break
        else:
            print(msg)


@common.auth('teacher')
def check_stu_from_course():
    while True:
        flag, course_list = teacher_interface.check_course_interface(
            teacher_info.get('user')
        )
        if not flag:
            print(course_list)
            break

        for index2, course_name in enumerate(course_list):
            print(f'编号：{index2},学校名：{course_name}')

        choice2 = input('请输入选择的课程编号：').strip()
        if not choice2.isdigit():
            print('输入有误')
            continue
        choice2 = int(choice2)

        if choice2 not in range(len(course_list)):
            print('输入编号有误！')
            continue

        course_name = course_list[choice2]

        flag2,student_list = teacher_interface.get_student_interface(
            course_name,teacher_info.get('user')
        )

        if flag2:
            print(student_list)
            break
        else:
            print(student_list)
            break






@common.auth('teacher')
def change_score_from_student():
    while True:
        flag, course_list = teacher_interface.check_course_interface(
            teacher_info.get('user')
        )
        if not flag:
            print(course_list)
            break

        for index2, course_name in enumerate(course_list):
            print(f'编号：{index2},学校名：{course_name}')

        choice2 = input('请输入选择的课程编号：').strip()
        if not choice2.isdigit():
            print('输入有误')
            continue
        choice2 = int(choice2)

        if choice2 not in range(len(course_list)):
            print('输入编号有误！')
            continue

        course_name = course_list[choice2]

        flag2,student_list = teacher_interface.get_student_interface(
            course_name,teacher_info.get('user')
        )

        if not flag2:
            print(student_list)
            break

        for index2,student_name in enumerate(student_list):
            print(f'编号：{index2}  学生名：{student_name}')

        choice3 = input('请输入学生编号：').strip()

        choice3 = int(choice3)

        if choice3 not in range(len(course_list)):
            print('输入编号有误！')
            continue

        student_name = student_list[choice3]
        score = input('请输入需要修改的成绩：').strip()
        if not score.isdigit():
            continue
        score = int(score)
        flag3,msg = teacher_interface.change_score_interface(
            course_name,student_name,score,teacher_info['user']
        )

        if flag3:
            print(msg)
            break



func_dict={
    '1':login,
    '2':check_course,
    '3':choose_course,
    '4':check_stu_from_course,
    '5':change_score_from_student
}




def teacher_view():
    while True:
        print('''
        老师
        1.登录
        2.查看教授课程
        3.选择教授课程
        4.查看课程下学生
        5.修改学生分数
        ''')
        choice = input('请输入功能编号：').strip()
        if choice == 'q':
            break

        if choice not in func_dict:
            print('输入有误，请重新输入')
            continue

        func_dict.get(choice)()