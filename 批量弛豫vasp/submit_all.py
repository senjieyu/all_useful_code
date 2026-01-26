#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pathlib import Path

def submit_all(jobs_dir: Path, start: int, end: int, dry_run: bool):
    jobs_dir = jobs_dir.resolve()
    if not jobs_dir.exists():
        raise SystemExit(f"[ERROR] jobs_dir not found: {jobs_dir}")

    os.chdir(jobs_dir)

    if end <= 0:
        # 自动探测最大编号
        mx = 0
        for p in Path(".").glob("bsub_*.lsf"):
            try:
                mx = max(mx, int(p.stem.split("_")[1]))
            except Exception:
                pass
        end = mx

    for i in range(start, end + 1):
        f = Path(f"bsub_{i}.lsf")
        if not f.exists():
            print("[WARN] missing", f.name)
            continue
        cmd = f"bsub<{f.name}"
        if dry_run:
            print("[DRY]", cmd)
        else:
            print("[SUBMIT]", cmd)
            os.system(cmd)

if __name__ == "__main__":
    # ============ 直接在这里改参数 ============
    JOBS_DIR = Path("jobs_script")
    START = 1
    END = 0      # 0 表示自动检测 bsub_*.lsf 的最大编号
    DRY_RUN = False
    # ========================================

    submit_all(JOBS_DIR, START, END, DRY_RUN)

