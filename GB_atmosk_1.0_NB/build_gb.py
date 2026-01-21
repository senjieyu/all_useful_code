#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
build_atomsk.py
功能：
- 自动创建工作目录 WORK_DIR 及每个 job 子目录
- 自动调用 Atomsk 生成 crystal1/crystal2 并 merge 成晶界
- 生成 polycrystal（写 polyX.txt 并调用 --polycrystal）
- 过程文件保留在 WORK_DIR/<jobname>/
- 最终结构统一复制到 FINAL_DIR/
- 关键修复：传给 Atomsk 的所有输入/输出文件都用绝对路径，避免 "file does not exist" 交互报错
"""

from __future__ import annotations

import importlib
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple


# ---------------------------
# 基础工具
# ---------------------------
def run_cmd(cmd: List[str], cwd: Path) -> None:
    """执行命令并打印，cwd 为工作目录。"""
    print("\n>>>", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=str(cwd))


def ensure_atomsk(atomsk: str) -> str:
    """
    检查 atomsk 是否可用：
    - 若给了绝对/相对路径（包含 / 或 \\），则检查文件存在
    - 否则检查 PATH 中是否能找到
    """
    if os.path.isabs(atomsk) or ("/" in atomsk) or ("\\" in atomsk):
        p = Path(atomsk)
        if not p.exists():
            raise FileNotFoundError(f"ATOMSK 不存在：{atomsk}")
        return str(p)
    if shutil.which(atomsk) is None:
        raise FileNotFoundError("找不到 atomsk：请加入 PATH 或在 settings.py 里填 ATOMSK 绝对路径")
    return atomsk


def auto_name(prefix: str, ext: str = ".cfg", **kv) -> str:
    """
    自动生成文件名，避免冲突、便于溯源。
    例如：Cu_tilt210001_dup=1x5x1.cfg
    """
    parts = [prefix]
    for k, v in kv.items():
        if v is None:
            continue
        if isinstance(v, tuple):
            v = "x".join(str(x).replace(".", "p") for x in v)
        else:
            v = str(v).replace(".", "p")
        parts.append(f"{k}={v}")
    return "_".join(parts) + ext


def export_final(work_out: Path, final_dir: Path, final_name: str) -> Path:
    """把工作目录里的最终输出复制到 FINAL_DIR。"""
    final_dir.mkdir(parents=True, exist_ok=True)
    dst = final_dir / final_name
    shutil.copy2(work_out, dst)
    print(f">>> Final exported -> {dst}")
    return dst


# ---------------------------
# Atomsk 封装：create / merge / poly
# （全部使用 .resolve() 绝对路径，彻底避免你遇到的报错）
# ---------------------------
def create_orient(atomsk: str, cwd: Path, lattice: str, a: float, element: str,
                 orient_xyz: List[str], dup: Tuple[int, int, int], out_file: Path) -> None:
    """
    对应官网：
      atomsk --create fcc a Cu orient X Y Z -duplicate nx ny nz crystal1.cfg
    """
    nx, ny, nz = dup
    cmd = [
        atomsk, "--create", lattice, f"{a}", element,
        "orient", orient_xyz[0], orient_xyz[1], orient_xyz[2],
        "-duplicate", str(nx), str(ny), str(nz),
        str(out_file.resolve()),  # ✅ 绝对路径
    ]
    run_cmd(cmd, cwd)


def create_rotate(atomsk: str, cwd: Path, lattice: str, a: float, element: str,
                 axis: str, angle: float, orth: bool,
                 dup: Tuple[int, int, int], out_file: Path) -> None:
    """
    对应官网：
      atomsk --create fcc a Cu -rotate Z 5.1 -orthogonal-cell -duplicate ... crystal1.cfg
    """
    nx, ny, nz = dup
    cmd = [
        atomsk, "--create", lattice, f"{a}", element,
        "-rotate", axis.upper(), f"{angle}",
    ]
    if orth:
        cmd.append("-orthogonal-cell")
    cmd += ["-duplicate", str(nx), str(ny), str(nz), str(out_file.resolve())]  # ✅
    run_cmd(cmd, cwd)


def merge_bicrystal(atomsk: str, cwd: Path, direction: str,
                    c1: Path, c2: Path, out_file: Path,
                    shift_sysid2: Tuple[float, float, float], wrap: bool) -> None:
    """
    对应官网：
      atomsk --merge Y 2 crystal1.cfg crystal2.cfg out.cfg
      可选：-select prop sysID 2 -shift dx dy dz
      可选：-wrap
    """
    dx, dy, dz = shift_sysid2
    cmd = [
        atomsk, "--merge", direction.upper(), "2",
        str(c1.resolve()),        # ✅ 绝对路径
        str(c2.resolve()),        # ✅
        str(out_file.resolve()),  # ✅
    ]
    if abs(dx) > 1e-12 or abs(dy) > 1e-12 or abs(dz) > 1e-12:
        cmd += ["-select", "prop", "sysID", "2", "-shift", f"{dx}", f"{dy}", f"{dz}"]
    if wrap:
        cmd.append("-wrap")
    run_cmd(cmd, cwd)


def write_poly(poly_path: Path, box: Tuple[float, float, float], rot_deg: float) -> None:
    """
    写 polycrystal 参数文件 polyX.txt（两晶粒示例）
    """
    Lx, Ly, Lz = box
    text = "\n".join([
        f"box {Lx} {Ly} {Lz}",
        f"node 0.5*box 0.25*box 0  0° 0° {-rot_deg}°",
        f"node 0.5*box 0.75*box 0  0° 0° {+rot_deg}°",
        ""
    ])
    poly_path.write_text(text, encoding="utf-8")
    print(f">>> poly file written: {poly_path}\n{text}")


def build_polycrystal(atomsk: str, cwd: Path, lattice: str, a: float, element: str,
                     box: Tuple[float, float, float], rot_deg: float,
                     poly_filename: str, out_file: Path) -> None:
    """
    自动：
    1) --create 生成基准晶胞 unitcell.xsf
    2) 写 polyX.txt
    3) --polycrystal 生成多晶
    """
    unitcell = cwd / "unitcell.xsf"
    polyfile = cwd / poly_filename

    run_cmd([atomsk, "--create", lattice, f"{a}", element, str(unitcell.resolve())], cwd)
    write_poly(polyfile, box, rot_deg)
    run_cmd([atomsk, "--polycrystal", str(unitcell.resolve()), str(polyfile.resolve()), str(out_file.resolve())], cwd)


# ---------------------------
# 主流程（全自动创建目录 + 执行开关）
# ---------------------------
def main() -> int:
    s = importlib.import_module("settings")

    atomsk = ensure_atomsk(s.ATOMSK)
    lattice = s.MATERIAL["lattice"]
    a = float(s.MATERIAL["a"])
    element = s.MATERIAL["element"]

    work_root = Path(s.WORK_DIR)
    final_dir = Path(s.FINAL_DIR)
    work_root.mkdir(parents=True, exist_ok=True)   # ✅ 自动创建
    final_dir.mkdir(parents=True, exist_ok=True)   # ✅ 自动创建

    wrap = bool(s.WRAP)
    shift2 = tuple(s.SHIFT_SYSID2)

    # 1) {210}[001] 对称倾转（orient）
    if bool(s.MAKE_TILT_210_001):
        job = "Cu_210_001_symmetric_tilt_by_orient"
        job_dir = work_root / job
        job_dir.mkdir(parents=True, exist_ok=True)  # ✅ 自动创建

        dup = tuple(s.DUPLICATE_TILT)
        c1 = job_dir / "crystal1.cfg"
        c2 = job_dir / "crystal2.cfg"
        out_work = job_dir / "bicrystal.cfg"

        # 官网示例取向
        o1 = ["[-120]", "[210]", "[001]"]
        o2 = ["[-120]", "[-2-10]", "[001]"]

        create_orient(atomsk, job_dir, lattice, a, element, o1, dup, c1)
        create_orient(atomsk, job_dir, lattice, a, element, o2, dup, c2)
        merge_bicrystal(atomsk, job_dir, "Y", c1, c2, out_work, shift2, wrap)

        final_name = auto_name(f"{element}_tilt210001_orient", dup=dup)
        export_final(out_work, final_dir, final_name)

    # 2) [001] 对称扭转（orient）
    if bool(s.MAKE_TWIST_001):
        job = "Cu_001_symmetric_twist_by_orient"
        job_dir = work_root / job
        job_dir.mkdir(parents=True, exist_ok=True)

        dup = tuple(s.DUPLICATE_TWIST)
        c1 = job_dir / "crystal1.cfg"
        c2 = job_dir / "crystal2.cfg"
        out_work = job_dir / "bicrystal.cfg"

        o1 = ["[-120]", "[001]", "[210]"]
        o2 = ["[-120]", "[001]", "[-2-10]"]

        create_orient(atomsk, job_dir, lattice, a, element, o1, dup, c1)
        create_orient(atomsk, job_dir, lattice, a, element, o2, dup, c2)
        merge_bicrystal(atomsk, job_dir, "Y", c1, c2, out_work, shift2, wrap)

        final_name = auto_name(f"{element}_twist001_orient", dup=dup)
        export_final(out_work, final_dir, final_name)

    # 3) 任意角对称倾转（rotate Z）
    if bool(s.MAKE_ARBITRARY_TILT):
        job = "Cu_arbitrary_symmetric_tilt_rotateZ"
        job_dir = work_root / job
        job_dir.mkdir(parents=True, exist_ok=True)

        dup = tuple(s.DUPLICATE_ARBITRARY_TILT)
        ang = float(s.ARBITRARY_ANGLE_DEG)

        c1 = job_dir / "crystal1.cfg"
        c2 = job_dir / "crystal2.cfg"
        out_work = job_dir / "bicrystal.cfg"

        create_rotate(atomsk, job_dir, lattice, a, element, "Z", +ang, True, dup, c1)
        create_rotate(atomsk, job_dir, lattice, a, element, "Z", -ang, True, dup, c2)
        merge_bicrystal(atomsk, job_dir, "Y", c1, c2, out_work, shift2, wrap)

        final_name = auto_name(f"{element}_tilt_rotZ", alpha=ang, dup=dup)
        export_final(out_work, final_dir, final_name)

    # 4) 任意角对称扭转（rotate Y）
    if bool(s.MAKE_ARBITRARY_TWIST):
        job = "Cu_arbitrary_symmetric_twist_rotateY"
        job_dir = work_root / job
        job_dir.mkdir(parents=True, exist_ok=True)

        dup = tuple(s.DUPLICATE_ARBITRARY_TWIST)
        ang = float(s.ARBITRARY_ANGLE_DEG)

        c1 = job_dir / "crystal1.cfg"
        c2 = job_dir / "crystal2.cfg"
        out_work = job_dir / "bicrystal.cfg"

        create_rotate(atomsk, job_dir, lattice, a, element, "Y", +ang, True, dup, c1)
        create_rotate(atomsk, job_dir, lattice, a, element, "Y", -ang, True, dup, c2)
        merge_bicrystal(atomsk, job_dir, "Y", c1, c2, out_work, shift2, wrap)

        final_name = auto_name(f"{element}_twist_rotY", alpha=ang, dup=dup)
        export_final(out_work, final_dir, final_name)

    # 5) polycrystal 两晶粒
    if bool(s.MAKE_POLYCRYSTAL_2G):
        job = "Cu_polycrystal_two_grains"
        job_dir = work_root / job
        job_dir.mkdir(parents=True, exist_ok=True)

        box = tuple(s.POLY_BOX)
        rot = float(s.POLY_ROT_DEG)
        out_work = job_dir / "polycrystal.cfg"

        build_polycrystal(atomsk, job_dir, lattice, a, element, box, rot, s.POLY_FILE_NAME, out_work)

        final_name = auto_name(f"{element}_poly2g", box=box, rot=rot)
        export_final(out_work, final_dir, final_name)

    print("\nAll done.")
    print(f"Final: {final_dir.resolve()}")
    print(f"Work:  {work_root.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
