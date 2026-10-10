# `os` 模块是 Python 中的标准库模块之一，它提供了与操作系统
# 交互的功能，可以让你执行文件和目录操作，管理进程，获取系统信
# 息等。以下是一些常用的 `os` 模块函数和方法：
#
# 1. **文件和目录操作**：
#
#    - `os.getcwd()`: 获取当前工作目录的路径。
#    - `os.chdir(path)`: 更改当前工作目录到指定的路径。
#    - `os.listdir(path)`: 返回指定目录中的文件和子目录的列表。
#    - `os.mkdir(path)`: 创建一个新目录。
#    - `os.makedirs(path)`: 递归创建多层目录。
#    - `os.remove(path)`: 删除指定文件。
#    - `os.rmdir(path)`: 删除指定目录。
#    - `os.removedirs(path)`: 递归删除多层目录。
#    - `os.rename(src, dst)`: 重命名文件或目录。
#
# 2. **路径操作**：
#
#    - `os.path.join(path, *paths)`: 将多个路径组合成一个。
#    - `os.path.exists(path)`: 判断路径是否存在。
#    - `os.path.isfile(path)`: 判断路径是否是一个文件。
#    - `os.path.isdir(path)`: 判断路径是否是一个目录。
#       os.path.getsize(path):查看文件的大小，字节
#       os.path.dirname()

#       os.path.basename()
#
# 3. **系统信息**：
#
#    - `os.name`: 返回操作系统的名称，例如 "posix" 或 "nt"。
#    - `os.uname()`: 返回关于当前系统的信息（仅在类 Unix 系统上可用）。
#
# 4. **进程管理**：
#
#    - `os.system(command)`: 在子 shell 中执行命令。
#    - `os.popen(command)`: 打开一个管道用于执行命令。
#    - `os.spawn*` 系列函数：用于创建新进程。
#    - `os.kill(pid, signal)`: 终止指定进程。
#
# 5. **环境变量**：
#
#    - `os.environ`: 包含系统环境变量的字典。
#
# 6. **文件操作**：
#
#    - `os.open(file, flags[, mode])`: 打开文件并返回文件描述符。
#    - `os.close(fd)`: 关闭文件描述符。
#    - `os.read(fd, n)`: 从文件描述符中读取数据。
#    - `os.write(fd, string)`: 向文件描述符写入数据。
#
# # 这只是 `os` 模块提供的一些常见功能。根据你的需求，你可以使用这些函
# 数来与操作系统进行各种交互，执行文件和目录操作，管理进程，获取系统信息等。