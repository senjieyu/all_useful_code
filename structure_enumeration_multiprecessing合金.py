import os
import time
import shutil
import itertools
from pathlib import Path
from functools import partial
from multiprocessing import Pool, cpu_count
from ase.io import read, write
from ase import Atoms
from icet.tools import enumerate_structures


def process_task(task, template_structure: Atoms, output_dir: str, tol: float = 1e-5):
    """
    单个任务：对某个二元体系 (A,B) + 固定 n_B 枚举所有不等价构型并写文件
    """
    pair_id, A, B, nB = task
    nsites = len(template_structure)

    if nB <= 0 or nB >= nsites:
        # A15B1...A1B15：nB只允许 1..nsites-1
        return (pair_id, f"{A}-{B}", nB, 0)

    concB = nB / nsites
    concentration_restrictions = {
        B: (concB - tol, concB + tol)
    }

    species = [A, B]

    # 输出目录：output_dir / "A-B" / "A{nA}B{nB}"
    nA = nsites - nB
    pair_dir = Path(output_dir) / f"{A}-{B}"
    comp_dir = pair_dir / f"{A}{nA}_{B}{nB}"
    comp_dir.mkdir(parents=True, exist_ok=True)

    print(f"[PID {os.getpid()}] 开始: ({A}-{B})  {A}{nA}{B}{nB}  conc({B})={concB:.5f}")

    start = time.time()
    enumerator = enumerate_structures(
        structure=template_structure,
        sizes=[1],
        chemical_symbols=species,
        concentration_restrictions=concentration_restrictions
    )

    count = 0
    for i, s in enumerate(enumerator, start=1):
        count += 1
        # 文件名：A-B__A{nA}B{nB}__id_0001.vasp
        filename = f"{A}-{B}__{A}{nA}{B}{nB}__id_{i:04d}.vasp"
        filepath = comp_dir / filename
        write(filepath, s, format="vasp", vasp5=True, direct=True, sort=True)

    elapsed = time.time() - start
    print(f"[PID {os.getpid()}] 完成: ({A}-{B}) {A}{nA}{B}{nB} -> {count} 个 (耗时 {elapsed:.2f}s)")

    return (pair_id, f"{A}-{B}", nB, count)


if __name__ == "__main__":
    # ----------------------------
    # 0) 你要枚举的元素集合（按图里的例子）
    # ----------------------------
    elements = ["Hf", "Nb", "Ti", "V", "Al"]

    # ----------------------------
    # 1) 读入结构并构造 BCC 2×2×2（你也可以换成你想要的）
    # ----------------------------
    atoms = read("POSCAR")
    template_structure = atoms * (2, 2, 2)   # 图里是 BCC 2*2*2 -> 常见为16位点
    nsites = len(template_structure)
    print(f"模板位点数: {nsites}")

    if nsites < 2:
        raise RuntimeError("位点数太少，无法做二元 A(n-1)B1 枚举。")

    # ----------------------------
    # 2) 输出目录
    # ----------------------------
    output_dir = "BinarySeeds_BCC_2x2x2"
    if os.path.exists(output_dir):
        print(f"清理旧目录: {output_dir}")
        shutil.rmtree(output_dir)
    Path(output_dir).mkdir(exist_ok=True)

    # ----------------------------
    # 3) 生成所有二元体系 + 成分空间 A(n-1)B1...A1B(n-1)
    # ----------------------------
    pairs = list(itertools.combinations(elements, 2))  # 不重复二元组合
    tasks = []
    pair_id = 0
    for (A, B) in pairs:
        pair_id += 1
        # 成分空间：B个数从 1 到 nsites-1
        for nB in range(1, nsites):
            tasks.append((pair_id, A, B, nB))

    print("\n--- 任务信息 ---")
    print(f"二元体系数量: {len(pairs)}")
    print(f"每个体系成分点数量: {nsites-1} (A{nsites-1}B1...A1B{nsites-1})")
    print(f"总任务数: {len(tasks)}")

    # ----------------------------
    # 4) 并行运行
    # ----------------------------
    num_processes = min(cpu_count(), len(tasks))
    print(f"使用进程池大小: {num_processes}\n")

    worker = partial(process_task, template_structure=template_structure, output_dir=output_dir)

    t0 = time.time()
    with Pool(processes=num_processes) as pool:
        results = pool.map(worker, tasks)
    t1 = time.time()

    # ----------------------------
    # 5) 汇总统计（输出一个 summary.txt）
    # ----------------------------
    total = sum(r[3] for r in results)
    summary_path = Path(output_dir) / "summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"Template sites = {nsites}\n")
        f.write(f"Elements = {elements}\n")
        f.write(f"Pairs = {pairs}\n\n")
        f.write("pair_id\tpair\tnB\tcount\n")
        for (pid, pair, nB, count) in results:
            f.write(f"{pid}\t{pair}\t{nB}\t{count}\n")
        f.write(f"\nTOTAL_STRUCTURES\t{total}\n")
        f.write(f"TOTAL_TIME_SEC\t{t1 - t0:.2f}\n")

    print("\n--- 全部完成 ---")
    print(f"总耗时: {t1 - t0:.2f} 秒")
    print(f"总生成结构数: {total}")
    print(f"输出目录: {output_dir}")
    print(f"汇总文件: {summary_path}")
