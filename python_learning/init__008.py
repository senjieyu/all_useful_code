# import sys
# sys.path.append(r'C:\Users\senji\Desktop\python_learning\learning\23.9.16包的尝试\foo\aaa')
# 内包不加这个也行，只要在导入时写foo.aaa.x就行了，如果加了上面的
# 则可以直接写     from .d import *


# from foo.aaa.d import *
# from foo.aaa.e import *
# from foo.aaa.f import *


#相对导入
from .d import *
from .e import *
from .f import *

#包内的导入建议使用相对导入
#因为如果模块名变的话，包内部的绝对导入都要改，出包用绝对导入