#模块查找优先级
#1.内存（内置模块）
#2.从硬盘找,按照sys.path中存放的文件顺序以此查找要导入的模块
import sys
print((sys.path))

# sys.modules查看已经加载到内存中的模块
#导入不在同一个文件夹下的文件，临时添加到环境变量中
#sys.pat h.append(r"文件所在文件夹的绝对路径")