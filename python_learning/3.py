# :爬取搜狗首页的页面数据
import requests
if __name__ == '__main__':

    # 兹step1:指定url
    url = 'https://www.sogou.com/'
    # 兹steg2:发起请求
    # 兹get方法会返回介响应对象
    response =requests.get(url=url)
    # 兹sten3:获取响应数据text返回的是宝符患形式的响应数据
    page_text =response.text
    print(page_text)
    # 共step4:持久化存馈
    with open( './s0gou.html', 'w',encoding=' utf-8 ' ) as fp:
        fp.write(page_text)
    print('爬取数据结束！！！')
