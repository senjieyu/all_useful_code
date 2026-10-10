import requests
import json
if __name__ == '__main__':
    # 指定URL
    post_url = 'https://fanyi.baidu.com/sug'
    # post请求参数处理（同get请求一致
    word = input('enter a word:')
    data ={
        'kw':word
    }
    # UA伪装：将对应的User-Agent封装到一个字典中
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.3'
    }
    # 请求发送
    response = requests.post(url=post_url,data=data,headers = headers)
    # 获取响应数据:json()方法返回的是obj （如果确认响应数据是json类型的，才可以使用json（））
    # Content-Type:application/json

    dict_obj=response.json()
    # 持久化存储
    fileName = word + '.json'
    fp = open(fileName, 'w', encoding=' utf-8')
    json.dump(dict_obj,fp=fp, ensure_ascii = False)
    print('over')




