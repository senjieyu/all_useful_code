
import requests, time, joblib, random, pathlib
from DrissionPage import WebPage

steel_items = joblib.load("steel_items_links.jb")

item_dicts_path = pathlib.Path("item_dicts")
if item_dicts_path.exists():
    item_dicts = joblib.load(item_dicts_path)
else:
    item_dicts = {}

# 加一个变体的list
variants_path = pathlib.Path("variants_path")
if variants_path.exists():
    variants = joblib.load(variants_path)
else:
    variants = []

base_url = "https://www.makeitfrom.com"
page = WebPage('s')
for link_index, link in enumerate(steel_items):
    
    if link_index >= len(item_dicts) + 10:
        break
    
    if link in item_dicts.keys():
        print("已存在，pass")
        continue
    
    item_dict = {}
    page.get(base_url+link)
    
    variant_div = page.eles("tag:a@@name=Intro")[0].nexts(".:split")
    if len(variant_div) > 0:
        item_dicts[link] = ''
        variants.extend([ i.attrs["href"] for i in variant_div[0].eles("tag:a")])
        print("变体，pass")
        sleep_time = random.randint(5, 15)
        print(f"睡眠{sleep_time}")
        time.sleep(sleep_time)
        print("睡眠结束")
        continue
    
    h1_title = page.eles("tag:a@@name=Intro")
    if len(h1_title) > 0:
        item_dict["title"] = h1_title[0].next().text
    
    young_modulus = page.eles("tag:p@text()^Elastic")
    if len(young_modulus) > 0: 
        young_modulus_p = young_modulus[0].next("tag:p").text.split("GPa")[0]
        item_dict[young_modulus[0].text] = young_modulus_p
    
    elongation = page.eles("tag:p@text()^Elongation")
    if len(elongation) > 0: 
        elongation_p = elongation[0].next("tag:p").text.split("GPa")[0]
        item_dict[elongation[0].text] = elongation_p
    
    fatigue_strength = page.eles("tag:p@text()^Fatigue")
    if len(fatigue_strength) > 0:
        fatigue_strength_p = fatigue_strength[0].next("tag:p").text.split("MPa")[0]
        item_dict[fatigue_strength[0].text] = fatigue_strength_p
    
    tensile_uts = page.eles("tag:p@text()$(UTS)")
    if len(tensile_uts) > 0:
        tensile_uts_p = tensile_uts[0].next("tag:p").text.split("MPa")[0]
        item_dict[tensile_uts[0].text] = tensile_uts_p
    
    tensile_proof = page.eles("tag:p@text()$(Proof)")
    if len(tensile_proof) > 0:
        tensile_proof_p = tensile_proof[0].next("tag:p").text.split("MPa")[0]
        item_dict[tensile_proof[0].text] = tensile_proof_p
    
    pren = page.eles("tag:p@text()^PREN")
    if len(pren) > 0:
        pren_p = pren[0].next("tag:p").text.split("MPa")[0]
        item_dict[pren[0].text] = pren_p
    
    resilience = page.eles("tag:p@text():Modulus of Resilience")
    if len(resilience) > 0:
        resilience_p = resilience[0].next("tag:p").text.split("kJ/m3")[0]
        item_dict[resilience[0].text] = resilience_p
    
    ele_dict = {}
    for index, ele in enumerate(page.eles("tag:tr")):
        tds = ele.eles("tag:td")
        ele_dict[tds[0].eles("tag:span")[1].text] = tds[2].text
    
    item_dict["elements"] = ele_dict
    
    item_dicts[link] = item_dict
    
    sleep_time = random.randint(5, 15)
    print(f"睡眠{sleep_time}")
    time.sleep(sleep_time)
    print("睡眠结束")
    
    if link_index % 5 == 0:
        joblib.dump(item_dicts, item_dicts_path)
        joblib.dump(variants, variants_path)


joblib.dump(item_dicts, item_dicts_path)
joblib.dump(variants, variants_path)
