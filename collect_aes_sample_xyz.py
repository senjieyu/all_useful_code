import os
from ase.io import iread,write
from tqdm import tqdm

def find_files_with_exact_name(root_dir, filename):
    """
    查找指定文件名的文件路径（精确匹配）

    Args:
        root_dir: 要搜索的根目录
        filename: 要查找的完整文件名（包括扩展名）

    Returns:
        list: 匹配的文件路径列表
    """
    matched_files = []

    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file == filename:
                full_path = os.path.join(root, file)
                matched_files.append(full_path)

    return matched_files

def has_lost_atoms_in_log(file_path):
    """
    检查单个log文件中是否包含"Lost atoms:"字符

    Args:
        file_path: log文件的完整路径

    Returns:
        bool: 如果包含"Lost atoms:"返回True，否则返回False
    """
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if "Lost atoms:" in line:
                    return True
        return False
    except Exception as e:
        print(f"无法读取文件 {file_path}: {e}")
        return False

# 使用示例
if __name__ == "__main__":
    import sys
    root_directory = sys.argv[1] #"1"  # 替换为你的目录路径
    search_keyword = "scf_filter.xyz"  # 替换为你要查找的关键词

    results = find_files_with_exact_name(root_directory, search_keyword)
    count = 0
    temp = []
    for file_path in tqdm(results):
        #print(file_path)
        data = list(iread(file_path))
        count += len(data)
        temp += data
    print(f'all_sample_stru:{count}')
    write('ase_sample.xyz',temp,format='extxyz')
