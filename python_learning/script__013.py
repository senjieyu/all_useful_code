# 如何处理异常
# 语法上的错误 SyntaxError

# 逻辑上的错误

# # 若错误的机制是不可预知的
# try:
#     子代码块# 有可能抛出异常的代码
# except 异常类型1 as e：
#     pass
# except 异常类型2 as d：
#     pass
# except 异常类型3 as f：
#     pass
# ...
# else:
#     如果被检测的字代码块没有异常发生，则会执行else的子代码
# finally:
#     无论被检测的子代码块有无异常发生，都会执行finally的子代码


#抛出异常同级别的后续代码不会运行


# except (异常类型2,异常类型2) as d：
#     print(d)
#
# except Exception as d：# 万能异常
#     print(d)










