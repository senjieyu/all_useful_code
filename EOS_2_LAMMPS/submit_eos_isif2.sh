#!/bin/bash
#BSUB -J eos_mlip_isif2
#BSUB -q 33
#BSUB -n 40
#BSUB -R "span[ptile=40]"
#BSUB -o %J.out
#BSUB -e %J.err

set -euo pipefail

DATA_DIR="eos_data"
INFILE="relax_isif2.in"
MLIP_INI="mlip.ini"
MTP="current.mtp"

# ✅ 用带 MLIP 的 LAMMPS
LMP_BIN="/work/phy-luoyf/apps/lammps_sus2/interface-lammps-mlip-2/lmp_intel_cpu_intelmpi"

# 每个结构用多少核（minimize 通常 1 就够）
RANKS_PER_TASK=1

NPROC="${LSB_DJOB_NUMPROC:-40}"
MAX_CONCURRENT=$(( NPROC / RANKS_PER_TASK ))
[ "${MAX_CONCURRENT}" -lt 1 ] && MAX_CONCURRENT=1

echo "[INFO] Host: $(hostname)"
echo "[INFO] NPROC: ${NPROC}"
echo "[INFO] LAMMPS: ${LMP_BIN}"
echo "[INFO] MAX_CONCURRENT: ${MAX_CONCURRENT}"
echo "[INFO] DATA_DIR: ${DATA_DIR}"
echo "[INFO] INFILE: ${INFILE}"
echo "[INFO] MLIP_INI: ${MLIP_INI}"
echo "[INFO] MTP: ${MTP}"

# ---- checks ----
[ -d "${DATA_DIR}" ] || { echo "[ERROR] ${DATA_DIR} not found"; exit 1; }
[ -f "${INFILE}" ]   || { echo "[ERROR] ${INFILE} not found"; exit 1; }
[ -f "${MLIP_INI}" ] || { echo "[ERROR] ${MLIP_INI} not found"; exit 1; }
[ -f "${MTP}" ]      || { echo "[ERROR] ${MTP} not found"; exit 1; }
[ -x "${LMP_BIN}" ]  || { echo "[ERROR] LAMMPS not executable: ${LMP_BIN}"; exit 1; }
command -v mpirun >/dev/null 2>&1 || { echo "[ERROR] mpirun not found"; exit 1; }

# Intel MPI on LSF
export I_MPI_HYDRA_BOOTSTRAP=lsf

mkdir -p eos_logs eos_results eos_relaxed eos_dump eos_xyz
: > eos_results/EV_raw.txt

shopt -s nullglob
FILES=( ${DATA_DIR}/data_V*.lmp )
[ ${#FILES[@]} -gt 0 ] || { echo "[ERROR] no ${DATA_DIR}/data_V*.lmp"; exit 1; }

# concurrency control
wait_for_slot () {
  while [ "$(jobs -rp | wc -l)" -ge "${MAX_CONCURRENT}" ]; do
    if wait -n 2>/dev/null; then :; else sleep 1; fi
  done
}

run_one () {
  local f="$1"
  local base tag vs oname outfile logfile dumpfile xyzfile

  base="$(basename "$f")"
  tag="${base%.lmp}"
  vs="$(echo "$tag" | sed -n 's/.*V\([0-9]\+\.[0-9]\+\).*/\1/p')"
  [ -z "${vs}" ] && vs="NA"

  oname="eos_relaxed/${tag}_relaxed.data"
  outfile="eos_results/${tag}.dat"
  logfile="eos_logs/${tag}.log"
  dumpfile="eos_dump/${tag}_final.dump"
  xyzfile="eos_xyz/${tag}_final.extxyz"

  echo "[RUN] ${base}  Vscale=${vs}"

  mpirun -bootstrap lsf -np "${RANKS_PER_TASK}" "${LMP_BIN}" \
    -in "${INFILE}" \
    -log "${logfile}" \
    -var fname "${f}" \
    -var oname "${oname}" \
    -var outfile "${outfile}" \
    -var mlip_ini "${MLIP_INI}" \
    -var dumpname "${dumpfile}"

  # --- convert LAMMPS custom dump -> ASE extended xyz ---
  # outfile columns: pe vol lx ly lz pxx pyy pzz pxy pxz pyz  (pressure in bar)
  # dumpfile has ITEM headers; we parse last snapshot and write extxyz:
  #   N
  #   Lattice="..." Properties=species:S:1:pos:R:3:forces:R:3 virial="..." energy=... stress="..." pbc="T T T"
  python - <<PY
import re
from pathlib import Path

tag = "${tag}"
vs = "${vs}"
dumpfile = Path("${dumpfile}")
outfile  = Path("${outfile}")
xyzfile  = Path("${xyzfile}")

if not dumpfile.exists() or dumpfile.stat().st_size == 0:
    raise SystemExit(f"[ERROR] missing dumpfile: {dumpfile}")
if not outfile.exists() or outfile.stat().st_size == 0:
    raise SystemExit(f"[ERROR] missing outfile: {outfile}")

# --- read scalars from outfile ---
vals = outfile.read_text().strip().split()
# pe vol lx ly lz pxx pyy pzz pxy pxz pyz
pe, vol = float(vals[0]), float(vals[1])
lx, ly, lz = map(float, vals[2:5])
pxx, pyy, pzz, pxy, pxz, pyz = map(float, vals[5:11])

# Convert pressure(bar) -> eV/Å^3
# 1 eV/Å^3 = 1.602176634e11 Pa ; 1 bar = 1e5 Pa
bar_to_eva3 = 1e5 / 1.602176634e11

# ASE "stress" 常用：张量(3x3) 以 eV/Å^3 表示，并且通常为“张力”为正（与压力相反）
# 这里取 stress = -pressure_tensor
sxx = -pxx * bar_to_eva3
syy = -pyy * bar_to_eva3
szz = -pzz * bar_to_eva3
sxy = -pxy * bar_to_eva3
sxz = -pxz * bar_to_eva3
syz = -pyz * bar_to_eva3

# virial (eV) 常见定义：virial = - stress * V  （因为 stress = -P）
# 所以 virial = P * V；这里用对称张量 3x3，输出9个数
vxx = -sxx * vol
vyy = -syy * vol
vzz = -szz * vol
vxy = -sxy * vol
vxz = -sxz * vol
vyz = -syz * vol

# --- parse last snapshot from LAMMPS dump custom ---
txt = dumpfile.read_text().splitlines()

# Find last "ITEM: TIMESTEP"
idxs = [i for i,l in enumerate(txt) if l.startswith("ITEM: TIMESTEP")]
if not idxs:
    raise SystemExit("[ERROR] no TIMESTEP in dump")
i0 = idxs[-1]

# Next lines format:
# i0: ITEM: TIMESTEP
# i0+1: timestep
# i0+2: ITEM: NUMBER OF ATOMS
# i0+3: N
# i0+4: ITEM: BOX BOUNDS ...
# i0+5..7: bounds
# i0+8: ITEM: ATOMS id type x y z fx fy fz
# i0+9..: atoms
N = int(txt[i0+3].strip())

# atom lines start:
atoms_start = i0 + 9
atoms = []
for line in txt[atoms_start:atoms_start+N]:
    parts = line.split()
    # id type x y z fx fy fz
    _id = int(parts[0]); typ = int(parts[1])
    x,y,z = map(float, parts[2:5])
    fx,fy,fz = map(float, parts[5:8])
    atoms.append((_id, typ, x,y,z, fx,fy,fz))
atoms.sort(key=lambda t: t[0])

# Map type->species (你这里是 Cu 单元素；如需多元素可扩展)
type_to_species = {1: "Cu"}
def sp(typ):
    return type_to_species.get(typ, f"X{typ}")

# Lattice string (orthogonal cell)
lattice = f'{lx} 0.0 0.0 0.0 {ly} 0.0 0.0 0.0 {lz}'

# stress and virial strings (9 comps)
stress9 = f'{sxx} {sxy} {sxz} {sxy} {syy} {syz} {sxz} {syz} {szz}'
virial9 = f'{vxx} {vxy} {vxz} {vxy} {vyy} {vyz} {vxz} {vyz} {vzz}'

# Extended XYZ comment line
comment = (
    f'Lattice="{lattice}" '
    f'Properties=species:S:1:pos:R:3:forces:R:3 '
    f'virial="{virial9}" '
    f'label=1 energy={pe} '
    f'stress="{stress9}" free_energy={pe} '
    f'pbc="T T T"'
)

# write extxyz
with xyzfile.open("w") as w:
    w.write(f"{N}\n")
    w.write(comment + "\n")
    for (_id, typ, x,y,z, fx,fy,fz) in atoms:
        w.write(f"{sp(typ):<2s}  {x: .8f}  {y: .8f}  {z: .8f}  {fx: .8f}  {fy: .8f}  {fz: .8f}\n")

print(f"[OK] wrote {xyzfile}")
PY
}

# run all points (parallel)
for f in "${FILES[@]}"; do
  wait_for_slot
  run_one "${f}" &
done
wait
echo "[INFO] All LAMMPS runs finished."

# gather results (still use first two columns: pe vol)
for f in "${FILES[@]}"; do
  base="$(basename "$f")"
  tag="${base%.lmp}"
  vs="$(echo "$tag" | sed -n 's/.*V\([0-9]\+\.[0-9]\+\).*/\1/p')"
  [ -z "${vs}" ] && vs="NA"

  outfile="eos_results/${tag}.dat"
  logfile="eos_logs/${tag}.log"

  if [ ! -s "${outfile}" ]; then
    echo "[WARN] missing ${outfile} (check ${logfile})"
    continue
  fi

  pe="$(awk '{print $1}' "${outfile}")"
  vol="$(awk '{print $2}' "${outfile}")"
  echo "${vs} ${vol} ${pe}" >> eos_results/EV_raw.txt
done

# plot
python - <<'PY'
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("eos_results/EV_raw.txt", sep=r"\s+", header=None,
                 names=["Vscale","Volume_A3","Energy_eV"])
df = df.sort_values("Volume_A3").reset_index(drop=True)
df["Eshift"] = df["Energy_eV"] - df["Energy_eV"].min()
df.to_csv("eos_results/eos_points.csv", index=False)

plt.figure()
plt.plot(df["Volume_A3"], df["Eshift"], marker="o")
plt.xlabel("Volume (A^3)")
plt.ylabel("E - Emin (eV)")
plt.title("EOS (MLIP/MTP, ISIF=2 fixed cell, relax atoms)")
plt.grid(True)
plt.tight_layout()
plt.savefig("eos_results/E_vs_V.png", dpi=200)
print(df)
print("\n[OK] eos_results/E_vs_V.png  eos_results/eos_points.csv")
PY

echo "[DONE]"

