#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import math
import shutil
from pathlib import Path
from typing import List, Tuple

try:
    from ase.io import read as ase_read
except Exception:
    ase_read = None


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


def kspacing_by_cell(poscar_path: Path, base: float, min_k: float, max_k: float, scale: float) -> float:
    """
    (可选) 根据晶胞大小自适应 KSPACING：
    - 晶胞越大 => KSPACING 越大（k网格更稀）
    - 晶胞越小 => KSPACING 越小（k网格更密）
    这里用平均晶格常数做一个简单缩放：k = base * (a_avg / scale)

    需要 ASE 支持读取 POSCAR。
    """
    if ase_read is None:
        return base

    try:
        atoms = ase_read(str(poscar_path), format="vasp")
        cell = atoms.cell.lengths()  # (a,b,c)
        a_avg = float(sum(cell) / 3.0)
        k = base * (a_avg / scale)
        if k < min_k:
            k = min_k
        if k > max_k:
            k = max_k
        return float(k)
    except Exception:
        return base


def write_incar_relax(incar_path: Path, ncore: int, kspacing: float, kgamma: bool):
    """
    全弛豫（晶胞+原子）：ISIF=3 + IBRION=2 + NSW>0
    KPOINTS 不生成，直接用 INCAR 的 KSPACING/KGAMMA
    """
    kg = ".TRUE." if kgamma else ".FALSE."
    text = f"""# ---------- Basic ----------
ENCUT  = 500
PREC   = Accurate
LREAL  = Auto
LASPH  = .TRUE.
ADDGRID= .TRUE.

LWAVE  = .FALSE.
LCHARG = .FALSE.

# ---------- Electronic ----------
ISMEAR = 0
SIGMA  = 0.05
NELM   = 200
NELMIN = 6
EDIFF  = 1E-7
GGA    = PE

# ---------- Relax (cell + ions) ----------
IBRION = 2
ISIF   = 3
NSW    = 200
EDIFFG = -0.01

# ---------- Symmetry ----------
ISYM   = 0

# ---------- k-mesh from INCAR ----------
KSPACING = {kspacing:.3f}
KGAMMA   = {kg}

# ---------- Parallel ----------
NCORE  = {ncore}
"""
    incar_path.write_text(text, encoding="utf-8")


def gen_potcar_by_vaspkit(sample_dir: Path):
    """
    只生成 POTCAR（不生成 KPOINTS）
    你原脚本用：(echo 103)|vaspkit
    """
    cwd = os.getcwd()
    os.chdir(sample_dir)
    os.system(r'(echo 103) | vaspkit > /dev/null 2>&1')
    os.chdir(cwd)


def build_lsf_head(dir_index: int, queue: str, cores: int, ptile: int, vasp_export_line: str):
    """
    按你给的模板：
    #BSUB -q 33
    #BSUB -n 40
    #BSUB -R "span[ptile=40]"
    """
    vasp_export = (vasp_export_line.strip() + "\n") if vasp_export_line.strip() else ""
    head = f"""#!/bin/bash
#BSUB -J a_job_{dir_index}
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
    incar_template_name: str,
    ncore_incar: int,

    use_adaptive_kspacing: bool,
    kspacing_base: float,
    kspacing_min: float,
    kspacing_max: float,
    kspacing_scale_ang: float,
    kgamma: bool,

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

    # 写 INCAR 模板（如果你用自适应 KSPACING，也会对每个样本单独写 INCAR；模板仅用于参考）
    incar_tpl = Path.cwd() / incar_template_name
    write_incar_relax(
        incar_path=incar_tpl,
        ncore=ncore_incar,
        kspacing=kspacing_base,
        kgamma=kgamma,
    )

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

            # 计算该样本的 KSPACING（固定 or 自适应）
            if use_adaptive_kspacing:
                ks = kspacing_by_cell(
                    poscar_path=sample_path / "POSCAR",
                    base=kspacing_base,
                    min_k=kspacing_min,
                    max_k=kspacing_max,
                    scale=kspacing_scale_ang,
                )
            else:
                ks = kspacing_base

            # 写该样本 INCAR（直接包含 KSPACING/KGAMMA）
            incar_path = sample_path / "INCAR"
            write_incar_relax(incar_path, ncore_incar, ks, kgamma)

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

    # 生成一个 start_calc.py（可选，按你原来的习惯）
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
    print(f"  POSCAR dirs : {filter_dir}")
    print(f"  LSF scripts : {jobs_dir}")
    print(f"  Submit help : {Path.cwd() / 'start_calc.py'}")
    print("  Note: KPOINTS not generated; using KSPACING/KGAMMA in INCAR.")


if __name__ == "__main__":
    # ===================== 直接在这里改参数 =====================

    # 输入：你的 *.vasp 文件目录
    IN_DIR = Path("xyz2poscar")

    # 输出目录
    FILTER_DIR = Path("filter")
    JOBS_DIR = Path("jobs_script")

    # dir_1..dir_N 的数量（相当于你原来的 calc_dir_num）
    CALC_DIR_NUM = 7

    # INCAR 里 NCORE
    NCORE_INCAR = 4

    # KSPACING / KGAMMA：不生成 KPOINTS
    # 方案1：固定 KSPACING
    USE_ADAPTIVE_KSPACING = False
    KSPACING_BASE = 0.15
    KGAMMA = True

    # 方案2：按晶胞大小自适应（需要 ASE 能读 POSCAR）
    # kspacing = base * (a_avg / scale_ang)
    # 常用：scale_ang 取 ~10A 作为“参照晶胞尺寸”
    KSPACING_MIN = 0.10
    KSPACING_MAX = 0.30
    KSPACING_SCALE_ANG = 10.0

    # LSF 参数（按你给的模板）
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
        incar_template_name="INCAR.relax.template",
        ncore_incar=NCORE_INCAR,

        use_adaptive_kspacing=USE_ADAPTIVE_KSPACING,
        kspacing_base=KSPACING_BASE,
        kspacing_min=KSPACING_MIN,
        kspacing_max=KSPACING_MAX,
        kspacing_scale_ang=KSPACING_SCALE_ANG,
        kgamma=KGAMMA,

        queue=BSUB_QUEUE,
        bsub_cores=BSUB_CORES,
        bsub_ptile=BSUB_PTILE,
        vasp_export_line=VASP_EXPORT_LINE,

        run_vaspkit_potcar=RUN_VASPKIT_POTCAR,
    )

