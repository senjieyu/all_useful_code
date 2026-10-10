def auth(role):

    from core import admin,student,teacher
    def login_auth(func):
        def inner(*args, **kwargs):
            if role == 'admin':
                if admin.admin_info['user']:
                    res = func(*args, **kwargs)  # res用于，万一func有返回值的情况
                    return res
                else:
                    admin.login()

            if role == 'student':
                if student.student_info['urer']:
                    res = func(*args, **kwargs)  # res用于，万一func有返回值的情况
                    return res
                else:
                    student.login()


            if role == 'teacher':
                if teacher.teacher_info['urer']:
                    res = func(*args, **kwargs)  # res用于，万一func有返回值的情况
                    return res
                else:
                    teacher.login()
            else:
                print('当前识图没有权限')
        return inner  # 第一步，返回内部的函数

    return login_auth
