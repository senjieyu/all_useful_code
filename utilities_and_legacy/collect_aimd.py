import os
import sys
from tqdm import tqdm
from ase.io import iread,read,write
import numpy as np

def find_files_with_keyword(root_dir, keyword):
    """
    查找包含特定关键词的文件

    参数:
        root_dir: 要搜索的根目录
        keyword: 要查找的关键词

    返回:
        包含关键词的文件路径列表
    """
    matching_files = []

    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if keyword in file:
                full_path = os.path.join(root, file)
                matching_files.append(full_path)

    return matching_files

def read_vasp_xml(file,out_file):
    for a in iread(file):
        s = a.get_stress()
        six2nine = np.array([s[0], s[5], s[4], s[5], s[1], s[3], s[4], s[3], s[2]])

        virial = six2nine * -1 * a.get_volume()
        a.info['virial'] = virial
        #a.info['path'] = os.path.basename(file)
        write(out_file,a,format='extxyz',append=True)
if __name__== '__main__':
    root_dir = sys.argv[1]
    keyword = 'vasprun.xml'
    out_file = 'AIMD.xyz'
    file_list = find_files_with_keyword(root_dir, keyword)
    for file in tqdm(file_list):
        read_vasp_xml(file, out_file)


