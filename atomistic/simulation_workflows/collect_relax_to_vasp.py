#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import csv
import shutil
from pathlib import Path
from typing import Optional, Dict, List

import numpy as np
from tqdm import tqdm

try:
    from ase.io import read as ase_read, write as ase_write
except Exception:
    ase_read = None
    ase_write = None


def has_ok(sample_dir: Path) -> bool:
    return (sample_dir / "__ok__").exists()


def choose_relaxed_structure(sample_dir: Path) -> Optional[Path]:
    """Prefer CONTCAR, fallback POSCAR."""
    contcar = sample_dir / "CONTCAR"
    poscar = sample_dir / "POSCAR"
    if contcar.exists() and contcar.stat().st_size > 0:
        return contcar
    if poscar.exists() and poscar.stat().st_size > 0:
        return poscar
    return None


def parse_nelm(incar_path: Path, default: int = 60) -> int:
    pattern = re.compile(r"^\s*NELM\s*=\s*(\d+)", re.IGNORECASE)
    if not incar_path.exists():
        return default
    with incar_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = pattern.search(line)
            if m:
                return int(m.group(1))
    return default


def oszicar_electronic_converged(oszicar_path: Path, nelm: int) -> Optional[bool]:
    """
    Heuristic:
      - if last DAV/RMM iteration index == NELM -> likely NOT converged -> False
      - else -> True
      - unreadable -> None
    """
    if not oszicar_path.exists():
        return None
    last_tokens = None
    try:
        with oszicar_path.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if "DAV:" in line or "RMM:" in line:
                    last_tokens = line.split()
        if not last_tokens:
            return None
        it = int(last_tokens[1])
        return (it != int(nelm))
    except Exception:
        return None


def max_force_from_vasprun(sample_dir: Path) -> Optional[float]:
    """Max |F| from vasprun.xml (eV/A). Needs ASE."""
    if ase_read is None:
        return None
    vxml = sample_dir / "vasprun.xml"
    if not vxml.exists():
        return None
    try:
        atoms = ase_read(str(vxml), format="vasp-xml")
        f = atoms.get_forces()
        return float(np.linalg.norm(f, axis=1).max())
    except Exception:
        return None


def export_structure(src: Path, dst: Path, mode: str = "copy"):
    """
    mode:
      - copy: copy file directly
      - vasp: read+write by ASE for normalization
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    if mode == "copy" or ase_read is None or ase_write is None:
        shutil.copy2(src, dst)
        return
    try:
        atoms = ase_read(str(src), format="vasp")
        ase_write(str(dst), atoms, format="vasp")
    except Exception:
        shutil.copy2(src, dst)


# -------- time parsing --------
TIME_LINE_RE = re.compile(r"job_(\d+)_(\d+)\s+total_runtime:([0-9.]+)\s*s")


def parse_time_txt(dir_path: Path) -> Dict[int, float]:
    """
    Parse filter/dir_i/time.txt -> {sample_index: runtime_seconds}
    Line format (from your lsf):
      job_i_item total_runtime:xxx s
    """
    tfile = dir_path / "time.txt"
    mp: Dict[int, float] = {}
    if not tfile.exists():
        return mp
    try:
        with tfile.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                m = TIME_LINE_RE.search(line)
                if not m:
                    continue
                sample_idx = int(m.group(2))
                runtime = float(m.group(3))
                mp[sample_idx] = runtime  # keep last record
    except Exception:
        return mp
    return mp


def time_by_markers(sample_dir: Path) -> Optional[float]:
    """Fallback runtime = mtime(__ok__) - mtime(__start__)."""
    st = sample_dir / "__start__"
    ok = sample_dir / "__ok__"
    if not (st.exists() and ok.exists()):
        return None
    try:
        dt = ok.stat().st_mtime - st.stat().st_mtime
        if dt < 0:
            return None
        return float(dt)
    except Exception:
        return None


def collect_relax_to_vasp(
    filter_dir: Path,
    output_dir: Path,
    incar_path_for_nelm: Path,
    use_ok_flag: bool = True,
    check_electronic_convergence: bool = True,
    force_threshold: Optional[float] = None,
    export_mode: str = "copy",
):
    filter_dir = filter_dir.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    out_success = output_dir / "success"
    out_failed_no_ok = output_dir / "failed" / "no_ok"
    out_failed_no_conv = output_dir / "failed" / "no_converge"
    out_failed_high_force = output_dir / "failed" / "high_force"
    out_failed_missing = output_dir / "failed" / "missing_file"

    for p in [out_success, out_failed_no_ok, out_failed_no_conv, out_failed_high_force, out_failed_missing]:
        p.mkdir(parents=True, exist_ok=True)

    nelm = parse_nelm(incar_path_for_nelm, default=60)

    # dirs: filter/dir_*
    dir_list = [d for d in filter_dir.iterdir() if d.is_dir() and d.name.startswith("dir_")]

    def dir_key(p: Path):
        try:
            return int(p.name.split("_")[1])
        except Exception:
            return p.name

    dir_list = sorted(dir_list, key=dir_key)

    # pre-read time maps
    time_maps: Dict[str, Dict[int, float]] = {d.name: parse_time_txt(d) for d in dir_list}

    # list samples
    all_samples: List[Path] = []
    for d in dir_list:
        subs = [s for s in d.iterdir() if s.is_dir()]
        subs_sorted = sorted(subs, key=lambda p: int(p.name) if p.name.isdigit() else p.name)
        all_samples.extend(subs_sorted)

    total = len(all_samples)

    # counters
    success_exported = 0
    success_ok_or_file = 0

    failed_no_ok = 0
    failed_no_converge = 0
    failed_high_force = 0
    failed_missing_file = 0

    failed_list: List[str] = []
    rows: List[Dict[str, str]] = []

    # time report records
    time_records: List[Dict[str, object]] = []
    total_runtime_s = 0.0
    runtime_count = 0

    for sample_dir in tqdm(all_samples, desc="Collecting"):
        dir_name = sample_dir.parent.name
        sample_name = sample_dir.name
        sample_idx = int(sample_name) if sample_name.isdigit() else None

        ok_flag = has_ok(sample_dir)
        has_any_file = any((sample_dir / fn).exists() for fn in ["CONTCAR", "vasprun.xml", "OUTCAR", "OSZICAR"])
        if ok_flag or has_any_file:
            success_ok_or_file += 1

        # runtime (time.txt first, then marker fallback)
        runtime_s = None
        time_source = "NA"

        if sample_idx is not None:
            runtime_s = time_maps.get(dir_name, {}).get(sample_idx)
            if runtime_s is not None:
                time_source = "time.txt"

        if runtime_s is None:
            runtime_s = time_by_markers(sample_dir)
            if runtime_s is not None:
                time_source = "marker(__start__/__ok__)"

        time_records.append({
            "dir": dir_name,
            "sample": sample_name,
            "path": str(sample_dir),
            "runtime_s": runtime_s,
            "source": time_source
        })

        if runtime_s is not None:
            total_runtime_s += float(runtime_s)
            runtime_count += 1

        # structure file
        src_struct = choose_relaxed_structure(sample_dir)

        electronic_conv = None
        max_force = None
        reason = ""

        # rule 1: ok flag
        if use_ok_flag and (not ok_flag):
            failed_no_ok += 1
            reason = "NO_OK_FLAG"
            if src_struct is not None:
                export_structure(src_struct, out_failed_no_ok / f"{dir_name}_{sample_name}.vasp", mode=export_mode)
            failed_list.append(f"{sample_dir}  [{reason}]")
            rows.append({
                "dir": dir_name,
                "sample": sample_name,
                "path": str(sample_dir),
                "ok": str(ok_flag),
                "electronic_converged": "" if electronic_conv is None else str(electronic_conv),
                "max_force": "" if max_force is None else f"{max_force:.6f}",
                "runtime_s": "" if runtime_s is None else f"{runtime_s:.3f}",
                "exported": "False",
                "export_path": "",
                "reason": reason
            })
            continue

        # rule 2: missing structure
        if src_struct is None:
            failed_missing_file += 1
            reason = "MISSING_CONTCAR_POSCAR"
            failed_list.append(f"{sample_dir}  [{reason}]")
            rows.append({
                "dir": dir_name,
                "sample": sample_name,
                "path": str(sample_dir),
                "ok": str(ok_flag),
                "electronic_converged": "" if electronic_conv is None else str(electronic_conv),
                "max_force": "" if max_force is None else f"{max_force:.6f}",
                "runtime_s": "" if runtime_s is None else f"{runtime_s:.3f}",
                "exported": "False",
                "export_path": "",
                "reason": reason
            })
            continue

        # rule 3: electronic convergence
        if check_electronic_convergence:
            electronic_conv = oszicar_electronic_converged(sample_dir / "OSZICAR", nelm)
            if electronic_conv is False:
                failed_no_converge += 1
                reason = f"NO_ELECTRONIC_CONV(NELM={nelm})"
                export_structure(src_struct, out_failed_no_conv / f"{dir_name}_{sample_name}.vasp", mode=export_mode)
                failed_list.append(f"{sample_dir}  [{reason}]")
                rows.append({
                    "dir": dir_name,
                    "sample": sample_name,
                    "path": str(sample_dir),
                    "ok": str(ok_flag),
                    "electronic_converged": str(electronic_conv),
                    "max_force": "" if max_force is None else f"{max_force:.6f}",
                    "runtime_s": "" if runtime_s is None else f"{runtime_s:.3f}",
                    "exported": "False",
                    "export_path": "",
                    "reason": reason
                })
                continue

        # rule 4: force threshold (optional)
        if force_threshold is not None:
            max_force = max_force_from_vasprun(sample_dir)
            if max_force is not None and max_force > float(force_threshold):
                failed_high_force += 1
                reason = f"HIGH_FORCE(maxF={max_force:.4f} > {force_threshold})"
                export_structure(src_struct, out_failed_high_force / f"{dir_name}_{sample_name}.vasp", mode=export_mode)
                failed_list.append(f"{sample_dir}  [{reason}]")
                rows.append({
                    "dir": dir_name,
                    "sample": sample_name,
                    "path": str(sample_dir),
                    "ok": str(ok_flag),
                    "electronic_converged": "" if electronic_conv is None else str(electronic_conv),
                    "max_force": f"{max_force:.6f}",
                    "runtime_s": "" if runtime_s is None else f"{runtime_s:.3f}",
                    "exported": "False",
                    "export_path": "",
                    "reason": reason
                })
                continue

        # success export
        out_path = out_success / f"{dir_name}_{sample_name}.vasp"
        export_structure(src_struct, out_path, mode=export_mode)
        success_exported += 1

        rows.append({
            "dir": dir_name,
            "sample": sample_name,
            "path": str(sample_dir),
            "ok": str(ok_flag),
            "electronic_converged": "" if electronic_conv is None else str(electronic_conv),
            "max_force": "" if max_force is None else f"{max_force:.6f}",
            "runtime_s": "" if runtime_s is None else f"{runtime_s:.3f}",
            "exported": "True",
            "export_path": str(out_path),
            "reason": "SUCCESS"
        })

    failed_total = total - success_exported

    # report.csv
    report_csv = output_dir / "report.csv"
    with report_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["dir", "sample", "path", "ok", "electronic_converged", "max_force", "runtime_s", "exported", "export_path", "reason"]
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    # failed_list
    (output_dir / "failed_list.txt").write_text("\n".join(failed_list) + "\n", encoding="utf-8")

    # convergence criteria record
    criteria = [
        "===== Convergence / Success Criteria Used in This Collection =====",
        "",
        "1) Success(export) criteria:",
        f"   - use_ok_flag = {use_ok_flag}",
        "   - If use_ok_flag=True: '__ok__' must exist in sample folder",
        "   - Structure file must exist: CONTCAR preferred; fallback POSCAR",
        "",
        "2) Electronic convergence heuristic:",
        f"   - check_electronic_convergence = {check_electronic_convergence}",
        f"   - NELM read from: {incar_path_for_nelm} (NELM={nelm})",
        "   - Rule: OSZICAR last DAV/RMM iteration index != NELM => converged",
        "   - If last iteration index == NELM => NOT converged (likely max steps reached)",
        "",
        "3) Ionic/relaxation convergence:",
        "   - This collector does NOT parse OUTCAR for 'reached required accuracy'.",
        "   - It relies on '__ok__' marker produced by your run scripts.",
        "",
        "4) Force threshold filter (optional):",
        f"   - force_threshold = {force_threshold}",
        "   - If set: max|F| from vasprun.xml (ASE) must satisfy max|F| <= threshold",
        "",
        "5) Runtime (per-sample):",
        "   - Primary: filter/dir_i/time.txt lines 'job_i_idx total_runtime:xxx s'",
        "   - Fallback: mtime(__ok__) - mtime(__start__)",
        "",
    ]
    (output_dir / "convergence_criteria.txt").write_text("\n".join(criteria) + "\n", encoding="utf-8")

    # summary
    summary = [
        "===== Relax Collection Summary =====",
        f"filter_dir                 : {filter_dir}",
        f"output_dir                 : {output_dir}",
        f"total_samples              : {total}",
        f"success_exported           : {success_exported}",
        f"success_ok_or_file         : {success_ok_or_file}",
        f"failed_total               : {failed_total}",
        "",
        "----- Failure breakdown -----",
        f"failed_no_ok_flag          : {failed_no_ok}",
        f"failed_missing_structure   : {failed_missing_file}",
        f"failed_no_electronic_conv  : {failed_no_converge} (NELM={nelm})",
        f"failed_high_force          : {failed_high_force} (threshold={force_threshold})",
        "",
        f"Outputs:",
        f"  - success structures: {out_success}",
        f"  - failed structures : {output_dir/'failed'}",
        f"  - report            : {report_csv}",
        f"  - summary           : {output_dir/'summary.txt'}",
        f"  - failed list       : {output_dir/'failed_list.txt'}",
        f"  - criteria          : {output_dir/'convergence_criteria.txt'}",
    ]
    (output_dir / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")

    # -------- time_report.txt (per-sample + total) --------
    time_report = output_dir / "time_report.txt"

    def fmt_time(sec: float) -> str:
        m = sec / 60.0
        h = sec / 3600.0
        return f"{sec:.3f} s | {m:.3f} min | {h:.3f} hr"

    lines = []
    lines.append("===== Time Report (Relax Calculations) =====")
    lines.append(f"filter_dir : {filter_dir}")
    lines.append(f"count_with_time : {runtime_count} / {total}")
    lines.append("")
    lines.append("per-sample runtime:")
    lines.append("FORMAT: dir  sample  runtime(s|min|hr)  source  path")
    lines.append("-" * 120)

    def sort_key(r):
        try:
            d = int(str(r["dir"]).split("_")[1])
        except Exception:
            d = str(r["dir"])
        try:
            s = int(str(r["sample"]))
        except Exception:
            s = str(r["sample"])
        return (d, s)

    for r in sorted(time_records, key=sort_key):
        rs = r["runtime_s"]
        if rs is None:
            rt = "NA"
        else:
            rt = fmt_time(float(rs))
        lines.append(f'{str(r["dir"]):8s} {str(r["sample"]):8s} {rt:30s} {str(r["source"]):22s} {str(r["path"])}')

    lines.append("-" * 120)
    lines.append("TOTAL runtime (sum of available samples):")
    lines.append(fmt_time(total_runtime_s))
    lines.append("")
    lines.append("Notes:")
    lines.append("- source=time.txt: parsed from filter/dir_i/time.txt written by your LSF script")
    lines.append("- source=marker(__start__/__ok__): fallback from file mtime difference")
    lines.append("- NA: no time information found for that sample")
    lines.append("")

    time_report.write_text("\n".join(lines), encoding="utf-8")

    print("\n".join(summary))
    print(f"\n[INFO] Time report written: {time_report}")


if __name__ == "__main__":
    # ============ 修改接口（直接改这里） ============
    PWD = Path(os.getcwd())
    FILTER_DIR = PWD / "filter"
    OUTPUT_DIR = PWD / "output"

    INCAR_PATH_FOR_NELM = PWD / "INCAR.relax.template"  # 改成你实际存在的文件
    USE_OK_FLAG = True
    CHECK_ELECTRONIC_CONVERGENCE = True

    FORCE_THRESHOLD = None      # 例如 0.05；None=不筛
    EXPORT_MODE = "copy"        # copy 或 vasp
    # ===============================================

    collect_relax_to_vasp(
        filter_dir=FILTER_DIR,
        output_dir=OUTPUT_DIR,
        incar_path_for_nelm=INCAR_PATH_FOR_NELM,
        use_ok_flag=USE_OK_FLAG,
        check_electronic_convergence=CHECK_ELECTRONIC_CONVERGENCE,
        force_threshold=FORCE_THRESHOLD,
        export_mode=EXPORT_MODE,
    )

