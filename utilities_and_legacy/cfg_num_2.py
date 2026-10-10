#!/usr/bin/env python
import sys

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} train.cfg")
        sys.exit(1)

    filename = sys.argv[1]

    stru_num = 0               # 结构数（BEGIN_CFG 块数）
    type_set = set()           # 全局出现过的 type
    type_counts = {}           # 各 type 的原子总数

    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    nlines = len(lines)

    while i < nlines:
        line = lines[i].strip()

        # 每遇到一个 BEGIN_CFG 就是一个结构
        if line.startswith('BEGIN_CFG'):
            stru_num += 1
            natoms = None
            i += 1
            continue

        # 找 Size -> 下一行是原子个数
        if line.startswith('Size'):
            # 下一行应该是原子数
            j = i + 1
            # 跳过空行
            while j < nlines and lines[j].strip() == '':
                j += 1
            if j < nlines:
                natoms = int(lines[j].split()[0])
            i = j + 1
            continue

        # 找到 AtomData 行后，接下来的 natoms 行是原子数据
        if line.startswith('AtomData:'):
            if natoms is None:
                raise RuntimeError("Found 'AtomData:' but natoms is not set (no 'Size' block before?).")

            # AtomData 这一行的下一行开始是原子坐标
            i += 1
            for k in range(natoms):
                if i + k >= nlines:
                    raise RuntimeError("Unexpected end of file when reading AtomData block.")
                parts = lines[i + k].split()
                if len(parts) < 2:
                    continue
                # 行格式: id  type  x  y  z  fx  fy  fz
                atype = parts[1]   # 这里是你文件里的 0 / 1 / 2 等
                type_set.add(atype)
                type_counts[atype] = type_counts.get(atype, 0) + 1

            # 跳过刚刚读完的 natoms 行
            i += natoms
            continue

        i += 1

    print(f"stru_num: {stru_num}")
    print(f"type_list: {sorted(type_set)}")
    print(f"num_type: {len(type_set)}")
    print("type_counts:")
    for t in sorted(type_counts, key=lambda x: int(x)):
        print(f"  type {t}: {type_counts[t]} atoms")

if __name__ == '__main__':
    main()

