'''程序的入口'''

#添加环境变量
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))


from core import src
#开始执行项目函数
if __name__=='__main__':
    src.run()

    # 先执行用户视图层






