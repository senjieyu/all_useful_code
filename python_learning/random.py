# import random
#
# random.random(): 生成一个0到1之间的随机浮点数。
#
# import random
# random_num = random.random()
# random.randint(a, b): 生成一个在闭区间 [a, b] 内的随机整数。
#
# import random
# random_int = random.randint(1, 10)
# random.uniform(a, b): 生成一个在区间 [a, b] 内的随机浮点数。
#
# import random
# random_float = random.uniform(1.0, 5.0)
# random.choice(seq): 从序列 seq 中随机选择一个元素。
#
# import random
# my_list = [1, 2, 3, 4, 5]
# random_element = random.choice(my_list)
# random.shuffle(seq): 随机打乱序列 seq 中的元素顺序（注意，这个操作会修改原始序列）。
#
# import random
# my_list = [1, 2, 3, 4, 5]
# random.shuffle(my_list)
# random.sample(population, k): 从总体 population 中随机选择 k 个不重复的元素，返回一个列表。
#
# import random
# my_list = [1, 2, 3, 4, 5]
# random_sample = random.sample(my_list, 3)



# random.
#应用：随机验证码
# print(chr(65))#A
# print(ord('a'))#97
import random
res=''
for i in range(6):
    shuzi=str(random.randint(0,9))
    zimu=chr(random.randint(65,90))
    s1=random.choice([zimu,shuzi])
    res+=s1
print(res)




