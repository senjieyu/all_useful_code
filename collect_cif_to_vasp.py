#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import hashlib
import os
import shutil
from pathlib import Path

def short_hash(text: str, n: int = 8) -> str:
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:n]

def safe_name(stem: str, rel_path: str, suffix: str) -> str:
    """
    生成不易冲突的文件名：<stem>__<hash><suffix>
    """
    h = short_hash(rel_path)
    return f"{stem}__{h}{suffix}"

def ensure_unique_path(dst: Path) -> Path:
    """
    如果 dst 已存在，追加序号保证唯一
    """
    if not dst.exists():
        return dst
    stem = dst.stem
    suffix = dst.suffix
    parent = dst.parent
    idx = 1
    while True:
        cand = parent / f"{stem}__{idx}{suffix}"
        if not cand.exists():
            return cand
        idx += 1

def try_convert_cif_to_vasp(cif_path: Path, vasp_path: Path) -> None:
    """
    尝试将 cif 转成 vasp(POSCAR)。
    优先 pymatgen，失败则尝试 ASE。
    """
    # 1) pymatgen
    try:
        from pymatgen.core import Structure
        from pymatgen.io.vasp import Poscar

        struct = Structure.from_file(str(cif_path))
        poscar = Poscar(struct)
        poscar.write_file(str(vasp_path))
        return
    except ImportError:
        pass
    except Exception as e:
        # pymatgen 装了但解析失败，继续尝试 ASE
        pymatgen_err = e
    else:
        pymatgen_err = None

    # 2) ASE
    try:
        from ase.io import read, write
        atoms = read(str(cif_path))
        # vasp5=True 输出 VASP5 POSCAR 格式
        write(str(vasp_path), atoms, format="vasp", vasp5=True, direct=True)
        return
    except ImportError:
        raise RuntimeError(
            "未检测到可用转换库。请安装 pymatgen 或 ase：\n"
            "  pip install pymatgen\n"
            "或\n"
            "  pip install ase\n"
        )
    except Exception as e:
        msg = f"ASE 转换失败: {e}"
        if 'pymatgen_err' in locals() and pymatgen_err is not None:
            msg = f"pymatgen 解析失败: {pymatgen_err}\n" + msg
        raise RuntimeError(msg)

def main():
    root = Path(__file__).resolve().parent
    out_cif = root / "all_cif"
    out_vasp = root / "all_vasp"
    out_cif.mkdir(parents=True, exist_ok=True)
    out_vasp.mkdir(parents=True, exist_ok=True)

    total_found = 0
    cif_copied = 0
    vasp_written = 0
    failed = 0

    # 遍历所有 cif
    for cif in root.rglob("*.cif"):
        # 跳过输出目录，避免循环
        try:
            cif.relative_to(out_cif)
            continue
        except ValueError:
            pass
        try:
            cif.relative_to(out_vasp)
            continue
        except ValueError:
            pass

        if not cif.is_file():
            continue

        total_found += 1

        rel = cif.relative_to(root).as_posix()
        stem = cif.stem

        # 1) 收集 cif：复制到 all_cif
        cif_dst = out_cif / safe_name(stem, rel, ".cif")
        cif_dst = ensure_unique_path(cif_dst)
        shutil.copy2(cif, cif_dst)
        cif_copied += 1

        # 2) 转换成 vasp：写到 all_vasp
        vasp_dst = out_vasp / safe_name(stem, rel, ".vasp")
        vasp_dst = ensure_unique_path(vasp_dst)

        try:
            try_convert_cif_to_vasp(cif, vasp_dst)
            vasp_written += 1
            print(f"[OK] {cif} -> {vasp_dst}")
        except Exception as e:
            failed += 1
            print(f"[FAIL] {cif}\n       {e}\n")

    print("\n===== Summary =====")
    print(f"Scan root : {root}")
    print(f"Found cif : {total_found}")
    print(f"Copied cif: {cif_copied} -> {out_cif}")
    print(f"Wrote vasp: {vasp_written} -> {out_vasp}")
    print(f"Failed    : {failed}")

if __name__ == "__main__":
    main()
