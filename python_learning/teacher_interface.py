
from db import models

def check_course_interface(teacher_name):
    # 获取当前老师对象
    teacher_obj = models.Teacher.select(teacher_name)

    # course_list = teacher_obj.course_list_from_tea
    course_list = teacher_obj.show_course()

    if not course_list:
        return False,'老师没有选择课程'

    return True,course_list



def add_course_interface(course_name,teacher_name):
    teacher_obj = models.Teacher.select(teacher_name)

    course_list = teacher_obj.course_list_from_tea
    if course_name in course_list:
        return False,'该课程已存在'

    teacher_obj.add_course(course_name)
    return True,'添加课程成功！'

def get_student_interface(course_name,teacher_name):
    teacher_obj = models.Teacher.select(teacher_name)
    student_list = teacher_obj.get_student(course_name)

    if not student_list:
        return False,'学生没有选择该课程'
    return True,student_list


def change_score_interface( course_name,student_name,
                            score,teacher_name):

    teacher_obj = models.Teacher.select(teacher_name)
    teacher_obj.change_score(course_name,student_name,
                            score)
    return True,'修改分数成功'






