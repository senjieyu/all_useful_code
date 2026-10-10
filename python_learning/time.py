import datetime
import time

# #time
# #时间的三种格式
# #1.时间戳：从1970年到现在的秒数
# import time
# print(time.time())
# #用于时间间隔的计算
#
# #2.按照某种格式现实的时间：2022-3-1 11.11.00
# print(time.strftime('%Y-%m-%d %H:%M:%S %p'))#p代表上午下午
# #用于展示时间
#
# #3.结构化时间
# res=time.localtime()
# print(res)
# print(res.tm_mday)
# print(res.tm_year)
# #用于获取时间的某一部分
#
# import datetime
# print(datetime.datetime.now())
# #2023-09-18 14:05:15.564030
# print(datetime.datetime.now()+datetime.timedelta(days=3))
# #时间的加减

# #时间模块需要掌握的操作
# #1.时间格式的转换
# # struct_time->时间戳
# import time
# s_time=time.localtime()
# print(s_time)
# print(time.mktime(s_time))#将结构化时间转化成时间戳
#
# #时间戳->struct_time
# tp_time=time.time()
# print(time.localtime(tp_time))
#
# #世界标准时间
# time.gmtime()

#struct_time ->格式化的字符串形式的时间
# s_time=time.localtime()
# print(time.strftime('%Y-%m-%d %H:%M:%S',s_time))

#将str格式转换成---> struct_time格式，格式需要对应上
# print(time.strptime('1988-03-03 11:11:11','%Y-%m-%d %H:%M:%S'))#格式对应上


# # 掌握：format string<---->timestamp
# # '1988-03-03 11:11:11'   +7
# # format string--->struct_timetimestamp
# struct_time=time.strptime('1988-03-03 11:11:11','%Y-%m-%d %H:%M:%S')
# res=time.mktime(struct_time)+7*86400
# print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(res)))

# print(time.asctime())
# #Mon Sep 18 16:08:14 2023


print(datetime.datetime.fromtimestamp(33333333))




