#列表生成式
# l=['ss_dsb','alk','lxx_dsb','sds_dsb','sdao_dsb','sdq']
# nel=[]
# for name in l:
#     nel.append(name) if name.endswith('dsb') else ...
# print(nel)
#三元表达式是if

#生成式是for + if，，，不加if代表全都加（if True）
#因为本身就在列表里所以不需要append
# l=['ss_dsb','alk','lxx_dsb','sds_dsb','sdao_dsb','sdq']
# nel=[name for name in l if name.endswith('dsb')]
# nel
# print(nel)

#练习：将所有小写变成大写，所有的名字去掉dsb
# l=['ss_dsb','alk','lxx_dsb','sds_dsb','sdao_dsb','sdq']
# nel=[name.replace('_dsb','').upper() for name in l ]
# print(nel)

#字典生成式
# keys=['sd','sds','frg','afbw']
# values=['as',1,2,9]
# res={key : 6  for key in keys }
# print(res)

# items=[('aaaa',55),('ssss',66),('ddd',999)]
# res={v[0]:v[1] for v in items if v[0]!='ddd'}
# print(res)
# res={k:v for k,v in items if k!='ddd'}
# print(res)

#集合生成式
# keys=['sd','sds','frg','afbw']
# 集合 = {v for v in keys }
# print(集合)

#没有元祖生成式，元祖不可改

#生成器表达式
# g=(i for i in range(10) if i >3)
# print(g)
# print(next(g))
# print(next(g))
# print(next(g))
# print(next(g))
# print(next(g))

# with open ('11.txt',mode='rt',encoding='utf-8') as f:
#     list1=[]
#     for line in f:
#         l=len (line)
#         list1.append(l)
#     s=sum(list1)
#     print(s)
#简化
# with open ('11.txt',mode='rt',encoding='utf-8') as f:
#     list1=[len(v) for v in f]#如果文件有很多行的话，那么列表将非常大
#     s=sum(list1)
#     print(s)
#优化
with open ('11.txt',mode='rt',encoding='utf-8') as f:
    #写成生成器生成式  小括号
    shengchengqi_1=(len(v) for v in f)#如果文件有很多行的话，那么列表将非常大
    s=sum(shengchengqi_1) #sum对可迭代对象求和，一个个取值
    print(s)