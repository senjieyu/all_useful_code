import numpy as np
from ase import Atom, Atoms
from ase.build import bulk
from icet.tools import enumerate_structures
from ase.io import write,read
from pathlib import Path
import os
import shutil
import time
from multiprocessing import Pool, cpu_count
from functools import partial

def create_mock_structure() -> Atoms:
    """
    创建 AlN 纤锌矿结构 (Wurtzite) 的 3x3x2 超胞。
    """
    template = bulk('AlN', 'wurtzite', a=3.11, c=4.98)
    # 打印原始模板的晶胞参数以进行核对
    lengths = template.get_cell_lengths_and_angles()[:3]
    # 创建 3x3x2 超胞，总共 72 个原子 (36 Al, 36 N)
    atoms = template * (2, 2, 2) 
    return atoms

def process_batch(pair_count, al_structure, n_structure, output_dir):
    """
    处理特定 Sc 对数的枚举任务的工作函数。
    将被并行调用。
    """
    nsites_al = len(al_structure)
    species = ["Al", "Sc", ]
    
    n_sc = pair_count
    n_al = nsites_al - n_sc
    
    # 计算浓度
    conc_sc = n_sc / nsites_al
    
    # icet 需要浓度范围
    tol = 1e-5
    concentration_restrictions = {
        'Sc': (conc_sc - tol, conc_sc + tol),
    }
    
    print(f"[进程 {os.getpid()}] 开始处理: {n_sc}对 Sc (Al剩余 {n_al} 个)...")
    
    # 执行枚举
    enumerator = enumerate_structures(
        structure=al_structure,
        sizes=[1], 
        chemical_symbols=species,
        concentration_restrictions=concentration_restrictions
    )
    
    local_count = 0
    start_time = time.time()
    
    for i, enumerated_al in enumerate(enumerator):
        local_count += 1
        
        # 重组结构
        final_structure = enumerated_al.copy()
        final_structure.extend(n_structure)
        
        # 写入文件
        # 注意：为了避免并行写入冲突，使用 batch 前缀和局部索引
        # 文件名格式: batch_X_id_YYYY_...
        filename = f"batch_{pair_count}_id_{local_count:04d}_Al{n_al}_Sc{n_sc}.vasp"
        filepath = Path(output_dir) / filename
        write(filepath, final_structure, format='vasp', vasp5=True, direct=True,sort=True)
    
    elapsed = time.time() - start_time
    print(f"[进程 {os.getpid()}] 完成: {n_sc}对 Sc -> 找到 {local_count} 个结构 (耗时 {elapsed:.2f}s)")
    
    return local_count

if __name__ == '__main__':
    # 1. 初始化结构与目录
#    full_structure = create_mock_structure()
    atoms = read("POSCAR")
    full_structure = atoms * ( 2,2,1 )
    
    output_dir = "AlScN_enumerated_structures"
    if os.path.exists(output_dir):
        print(f"清理旧目录: {output_dir}")
        shutil.rmtree(output_dir)
    Path(output_dir).mkdir(exist_ok=True)

    # 2. 准备数据
    # 提取 Al 子晶格
    al_indices = [atom.index for atom in full_structure if atom.symbol == 'Al']
    al_structure = full_structure.copy()
    al_structure = al_structure[al_indices]
    
    # 提取 N 子晶格
    n_indices = [atom.index for atom in full_structure if atom.symbol == 'N']
    n_structure = full_structure.copy()
    n_structure = n_structure[n_indices]

    # 3. 设置并行任务
    max_pairs = 16
    # 创建任务列表：range(1, 6) 对应 1对 到 5对
    tasks = range(1, max_pairs + 1)
    
    # 获取可用 CPU 核心数 (根据任务量，最多用 max_pairs 个核心即可)
    num_processes = min(cpu_count(), len(tasks))
    print(f"\n--- 启动并行处理 ---")
    print(f"总任务数: {len(tasks)} (Sc 对数 1 到 {max_pairs})")
    print(f"使用进程池大小: {num_processes}")
    
    # 使用 partial 固定部分参数，只留 pair_count 变化
    worker_func = partial(process_batch, 
                          al_structure=al_structure, 
                          n_structure=n_structure, 
                          output_dir=output_dir)
    
    start_total = time.time()
    
    # 4. 执行并行池
    with Pool(processes=num_processes) as pool:
        # map 会阻塞直到所有结果返回，结果顺序与 tasks 一致
        results = pool.map(worker_func, tasks)
        
    end_total = time.time()
    
    # 5. 汇总结果
    total_structures = sum(results)
    print("\n--- 全部完成 ---")
    print(f"总耗时: {end_total - start_total:.2f} 秒")
    print(f"总共生成的结构数: {total_structures}")
    print(f"文件已保存在: {output_dir}")
