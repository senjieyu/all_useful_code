'''数据处理层'''
import os
import json
from conf.settings import BASE_PATH

def select(username):
    use_path=os.path.join(BASE_PATH,f'{username}.json')
    # 查看用户是否存在
    if os.path.exists(use_path):
        with open(use_path, 'r', encoding='utf-8') as f:
            user_dic = json.load(f)
            return user_dic


def save(user_dic):
    use_path = os.path.join(BASE_PATH, f'{user_dic["username"]}.json')
    # 写入用户
    with open(use_path, 'w', encoding='utf-8') as f:
        json.dump(user_dic, f, ensure_ascii=False)