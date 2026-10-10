
import requests, time, joblib, random
from bs4 import BeautifulSoup

steel_catos = joblib.load("steel_catos.jb")

base_url = "https://www.makeitfrom.com"
group_items = []
for index, i in enumerate(steel_catos):
    response = requests.get(base_url+i.attrs["href"])
    bs = BeautifulSoup(response.text.encode(), features="html.parser")
    group_item = bs.select("a[href*='material-properties']")
    links = [i.attrs["href"] for i in group_item]
    group_items.extend(links)
    time.sleep(random.randint(1, 8))
    # if index == 1:
    # break
    # print(index)

joblib.dump(group_items, "steel_items_links.jb")
