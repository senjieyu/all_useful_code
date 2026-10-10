#包就是一个包含有__init__.py文件的文件夹
#为何要有包
#包的本质是模块的一种形式，包是用来被当做模块导入
import sys
sys.path.append(r"C:\Users\senji\Desktop\python_learning\learning\23.9.15\foo")


# # from a1 import f1
# # from a2 import f2
# # from a3 import f3
# foo__init__文件中配合上述代码
# import foo
# foo.f2()


# # import a1
# # import a2
# # import a1
# import foo
# foo.a2.f2()

#.的左侧必须是个包，不能是文件

#相对导入
#.:代表当前文件夹
#..：代表上一级文件夹
#...:代表上上层文件夹
#以此类推







