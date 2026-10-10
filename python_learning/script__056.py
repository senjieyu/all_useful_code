#打补丁在首次导入的时候打
import json
import ujson
def monkey_patch_json():
    json.__name__ = 'ujson'
    json.dumps = ujson.dumps
    json.loads = ujson.loads
monkey_patch_json()

#第一次导入会导入到内存中
#后续再导入模块时，都是从内存中导入