#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import hashlib
import os
import shutil
from pathlib import Path

def short_hash(text: str, n: int = 8) -> str:
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:n]

def collect_vasp_files(root: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    count_total = 0
    count_copied = 0
    count_skipped_same = 0

    # 遍历所有 .vasp（递归）
    for src in root.rglob("*.vasp"):
        # 跳过输出目录自身，避免把复制过去的文件又扫描到
        try:
            src.relative_to(out_dir)
            continue
        except ValueError:
            pass

        if not src.is_file():
            continue

        count_total += 1

        rel = src.relative_to(root).as_posix()
        base = src.name  # 原文件名

        dst = out_dir / base

        # 如果目标已存在，避免覆盖：用“相对路径哈希”生成新名字
        if dst.exists():
            stem = src.stem
            suffix = src.suffix  # .vasp
            h = short_hash(rel)
            dst = out_dir / f"{stem}__{h}{suffix}"

            # 如果依然冲突，再追加序号
            idx = 1
            while dst.exists():
                dst = out_dir / f"{stem}__{h}__{idx}{suffix}"
                idx += 1

        # 若源和目标其实是同一文件（极少数情况），跳过
        try:
            if src.resolve() == dst.resolve():
                count_skipped_same += 1
                continue
        except Exception:
            pass

        shutil.copy2(src, dst)
        count_copied += 1
        print(f"[COPY] {src} -> {dst}")

    print("\n===== Done =====")
    print(f"Scan root: {root}")
    print(f"Output dir: {out_dir}")
    print(f"Found .vasp: {count_total}")
    print(f"Copied: {count_copied}")
    print(f"Skipped (same file): {count_skipped_same}")

def main():
    script_dir = Path(__file__).resolve().parent
    out_dir = script_dir / "all_sample"
    collect_vasp_files(script_dir, out_dir)

if __name__ == "__main__":
    main()
