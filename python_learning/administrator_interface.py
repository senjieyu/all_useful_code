from lib import common
from db import db_handler
import os

def renew():
    ad_dic={'school':[],
            'class':{},
            'teacher': {},
            'course':[],
            "username":'系统统一'}
    db_handler.save(ad_dic,'管理员')
if not os.path.exists(r'C:\Users\senji\Desktop\python_learning\learning\23.9.27项目作业\db\数据\系统统一管理员.pkl'):
    renew()

def login_interface(username,password):
        if  password=='123' and username=='asd':
            return True,'登录成功'
        else:
            return False,'密码错误'

# 学校类
class AdSet:
    ad_dic = db_handler.select('系统统一', '管理员')
    school = ad_dic['school']
    dict_class_name = ad_dic['class']
    def add_school(self,school_name):
        self.school.append(school_name)
        number = self.school.index(school_name)
        self.dict_class_name[self.school[number]] = {}
        db_handler.save(self.ad_dic,'管理员')
    def remove_school(self,school_name):
        ad_dic = db_handler.select('系统统一', '管理员')
        school = ad_dic['school']
        dict_class_name = ad_dic['class']
        school.pop(school_name)
        ad_dic['school'] = school
        number = school.index(school_name)
        del dict_class_name[school[number]]
        db_handler.save(ad_dic,'管理员')
    def show_school(self):
        ad_dic = db_handler.select('系统统一', '管理员')
        school = ad_dic['school']
        print(school)
        for s in school:
            print(s)


class CreatedClass:

    def chose_school (self,school_name):
        ad_dic = db_handler.select('系统统一', '管理员')
        school = ad_dic['school']
        if school_name in school:
            return True,f'您正在为{school_name}创建班级'
        else :
            return False,'请输入正确的学校名'

    def creat_class(self,school_name,class_name):
        ad_dic = db_handler.select('系统统一', '管理员')
        dict_class_name = ad_dic['class']   #{学校名：[班级，[老师]，[课程]]]   }
        list1=dict_class_name[school_name] = [1,2,3]   #班级里是列表   0为教师    1-为课程
        list1[0] = class_name
        list1[1] = []#老师
        list1[2] = []#课程
        db_handler.save(ad_dic, '管理员')
        res = input('继续添加请按 1 ，取消请按任意键：')
        if res == '1':
            return True,'您选择继续添加'
        else:
            return False,'您选择退出添加'

    def show_class_and_school(self):
        ad_dic = db_handler.select('系统统一', '管理员')
        school = ad_dic['school']
        dict_class_name = ad_dic['class']
        print(school)
        for s in school:
            print(s)
        print(dict_class_name)
        for a in dict_class_name.items():
            s,v =a[0],a[1]
            print(f'{s}:{v}')

    def del_class(self,school_name,class_name):
        ad_dic = db_handler.select('系统统一', '管理员')
        dict_class_name = ad_dic['class']
        if school_name in dict_class_name:
            dict_class_name[school_name].remove(class_name)
            db_handler.save(ad_dic, '管理员')
            return True,'操作成功'
        else:
            return False,'请输入正确的班级名'















