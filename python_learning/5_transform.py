
import requests, time, joblib, random, pathlib, pandas as pd
from collections import OrderedDict

item_dicts_path = pathlib.Path("item_dicts")
item_dicts = joblib.load(item_dicts_path)
item_dicts_path2 = pathlib.Path("item_dicts_for_variants.jb")
item_dicts2 = joblib.load(item_dicts_path2)
item_dicts.update(item_dicts2)

new_item_dicts = {}
for key, item in item_dicts.items():
    if item != '':
        new_item = OrderedDict()
        for k, i in item.items():
            if k != "elements":
                new_item[k] = i
        for a, b in item["elements"].items():
            new_item[a] = b
        new_item_dicts[key] = new_item

a = pd.DataFrame.from_dict(new_item_dicts).T
a.to_excel("makeitfrom_steels.xlsx")
