#获取所有学校名称接口
from conf import  settings
import os
from db import models
#获取所有学校名称接口
def get_all_school_interface():
    school_dir = os.path.join(
        settings.DB_PATH,'School'
    )
    if not os.path.exists(school_dir):
        return False,'没有学校，请先联系管理员'

    school_list = os.listdir(school_dir)
    return True,school_list


# 公共登录接口
def login_interface(user,pwd,user_type):

    if user_type == 'admin':
        obj = models.Admin.select(user)
    elif user_type == 'student':
        obj = models.Student.select(user)
    elif user_type == 'teacher':
        obj = models.Teacher.select(user)
    else:
        return False,'登录角色有误！'

    if obj:
        if pwd == obj.pwd :
            return True,'登录成功'
        else:
            return False,'密码错误'
    else:
        return False,'用户不存在！'



# 获取学校下所有课程接口
def get_course_in_school_interface(school_name):

    school_obj = models.School.select(school_name)

    course_list = school_obj.course_list

    if not course_list:
        return False,'该学校没有课程'
    return True,course_list
















