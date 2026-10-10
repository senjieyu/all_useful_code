import pickle
import os
from conf.setting import BASE_PATH
def select(username,role):
    use_path=os.path.join(BASE_PATH,'db','数据',f'{username}{role}.pkl')
    # 查看用户是否存在
    if os.path.exists(use_path):
        with open(use_path, 'rb') as f:
            user_dic = pickle.load(f)
            return user_dic

def save(user_dic,role):
    use_path = os.path.join(BASE_PATH,'db','数据',f'{user_dic["username"]}{role}.pkl')
    # 写入用户
    with open(use_path, 'wb') as f:
        pickle.dump(user_dic,f)