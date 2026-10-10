# In CPython, the global interpreter lock, or GIL, is a mutex that
# prevents multiple
# native threads from executing Python bytecodes at once. This lock is
# necessary mainly
# because CPython's memory management is not thread-safe. (However, since
# the GIL
# exists, other features have grown to depend on the guarantees that it
# enforces.)

# 重点:
# 1.GIL不是python的特点而是cPython解释器的特点
# 2.GIL是保证解释器级别的数据的安全
# 3.GIL会导致同一个进程下的多个线程的无法同时执行即无法利用多核优势(******)
# 4.针对不同的数据还是需要加不同的锁处理
# 5.解释型语言的通病:同一个进程下多个线程无法利用多核优势
#







