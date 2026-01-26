#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
from pathlib import Path
from typing import List, Tuple


def mkdir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def natural_sort_key(p: Path):
    # 1.vasp, 2.vasp, 10.vasp ...
    try:
        return int(p.stem)
    except Exception:
        return p.name


def split_ranges(num_poscar: int, calc_dir_num: int) -> List[Tuple[int, int]]:
    """把总样本数 num_poscar 均分到 calc_dir_num 个 dir 中，返回每个 dir 对应的 [start,end] (1-based)"""
    quotient = num_poscar // calc_dir_num
    remainder = num_poscar % calc_dir_num
    temp = [quotient] * calc_dir_num
    divide = [a + 1 if idx < remainder else a for idx, a in enumerate(temp)]
    ranges = [[sum(divide[:a + 1]) - divide[a] + 1, sum(divide[:a + 1])] for a in range(len(divide))]
    return [(r[0], r[1]) for r in ranges]


def gen_potcar_by_vaspkit(sample_dir: Path):
    """
    只生成 POTCAR（不生成 KPOINTS）
    原脚本用：(echo 103)|vaspkit
    """
    cwd = os.getcwd()
    os.chdir(sample_dir)
    os.system(r'(echo 103) | vaspkit > /dev/null 2>&1')
    os.chdir(cwd)


def build_lsf_head(dir_index: int, queue: str, cores: int, ptile: int, vasp_export_line: str):
    """
    按模板：
    #BSUB -q 33
    #BSUB -n 40
    #BSUB -R "span[ptile=40]"
    """
    vasp_export = (vasp_export_line.strip() + "\n") if vasp_export_line.strip() else ""
    head = f"""#!/bin/bash
#BSUB -J batch_relax_{dir_index}
#BSUB -q {queue}
#BSUB -n {cores}
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile={ptile}]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
{vasp_export}COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP vasp_std"
"""
    return head


def build_lsf_body(dir_index: int, range1: int, range2: int):
    body = f"""dir_1="filter"
dir_2="dir_{dir_index}"
cd ../$dir_1/$dir_2

path=$(pwd)

for item in {{{range1}..{range2}}}; do
    cd $path/$item
    start_time=$(date +%s.%N)
    touch __start__
    $COMMAND_std > vasp.out 2>&1
    touch __ok__
    end_time=$(date +%s.%N)
    runtime=$(echo "$end_time - $start_time" | bc)
    cd $path
    cd ..
    echo "job_{dir_index}_$item total_runtime:$runtime s" >> time.txt
done
"""
    return body


def build_relax_dirs(
    in_dir: Path,
    filter_dir: Path,
    jobs_dir: Path,
    calc_dir_num: int,

    queue: str,
    bsub_cores: int,
    bsub_ptile: int,
    vasp_export_line: str,

    run_vaspkit_potcar: bool,
):
    mkdir(filter_dir)
    mkdir(jobs_dir)

    # 收集 *.vasp
    files = sorted(in_dir.glob("*.vasp"), key=natural_sort_key)
    if not files:
        raise SystemExit(f"[ERROR] No *.vasp found in {in_dir}")

    num_poscar = len(files)
    calc_dir_num = min(calc_dir_num, num_poscar)
    ranges = split_ranges(num_poscar, calc_dir_num)

    # 统一 INCAR：使用“脚本所在目录”的 INCAR 作为模板
    script_dir = Path(__file__).resolve().parent
    incar_src = script_dir / "INCAR"
    if not incar_src.exists():
        raise SystemExit(f"[ERROR] INCAR not found next to script: {incar_src}")

    # 逐个样本建目录
    for dir_i in range(1, calc_dir_num + 1):
        dir_path = filter_dir / f"dir_{dir_i}"
        mkdir(dir_path)

        r1, r2 = ranges[dir_i - 1]
        for idx in range(r1, r2 + 1):
            sample_path = dir_path / f"{idx}"
            mkdir(sample_path)

            # POSCAR
            shutil.copy2(files[idx - 1], sample_path / "POSCAR")

            # INCAR：直接拷贝模板（不修改、不写入）
            shutil.copy2(incar_src, sample_path / "INCAR")

            # POTCAR by vaspkit
            if run_vaspkit_potcar:
                gen_potcar_by_vaspkit(sample_path)

    # 生成每个 dir 的 bsub 脚本
    for dir_i in range(1, calc_dir_num + 1):
        r1, r2 = ranges[dir_i - 1]
        lsf_path = jobs_dir / f"bsub_{dir_i}.lsf"
        head = build_lsf_head(dir_i, queue, bsub_cores, bsub_ptile, vasp_export_line)
        body = build_lsf_body(dir_i, r1, r2)
        lsf_path.write_text(head + "\n" + body, encoding="utf-8")

    # 生成 start_calc.py（按原习惯）
    start_calc = f"""import os

start = 1
end = {calc_dir_num}

pwd = os.getcwd()
jobs_script = os.path.join(pwd, "{jobs_dir.name}")
os.chdir(jobs_script)

for i in range(start, end+1):
    if not os.path.exists(f"bsub_{{i}}.lsf"):
        print("missing:", f"bsub_{{i}}.lsf")
        continue
    os.system(f"bsub<bsub_{{i}}.lsf")
"""
    (Path.cwd() / "start_calc.py").write_text(start_calc, encoding="utf-8")

    print("[OK] Done.")
    print(f"  INCAR template : {incar_src}")
    print(f"  POSCAR dirs    : {filter_dir}")
    print(f"  LSF scripts    : {jobs_dir}")
    print(f"  Submit help    : {Path.cwd() / 'start_calc.py'}")
    print("  Note: INCAR is copied as-is; no adaptive changes, no writing.")


if __name__ == "__main__":
    # ===================== 直接在这里改参数 =====================

    # 输入：你的 *.vasp 文件目录
    IN_DIR = Path("xyz2poscar")

    # 输出目录
    FILTER_DIR = Path("filter")
    JOBS_DIR = Path("jobs_script")

    # dir_1..dir_N 的数量（相当于你原来的 calc_dir_num）
    CALC_DIR_NUM = 7

    # LSF 参数
    BSUB_QUEUE = "33"
    BSUB_CORES = 40
    BSUB_PTILE = 40

    # 若集群需要指定 VASP 路径/环境（不需要就写 ""）
    VASP_EXPORT_LINE = "export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH"

    # 是否用 vaspkit 生成 POTCAR（需要你配置好赝势路径）
    RUN_VASPKIT_POTCAR = True

    # ==========================================================

    build_relax_dirs(
        in_dir=IN_DIR.resolve(),
        filter_dir=FILTER_DIR.resolve(),
        jobs_dir=JOBS_DIR.resolve(),
        calc_dir_num=CALC_DIR_NUM,

        queue=BSUB_QUEUE,
        bsub_cores=BSUB_CORES,
        bsub_ptile=BSUB_PTILE,
        vasp_export_line=VASP_EXPORT_LINE,

        run_vaspkit_potcar=RUN_VASPKIT_POTCAR,
    )
