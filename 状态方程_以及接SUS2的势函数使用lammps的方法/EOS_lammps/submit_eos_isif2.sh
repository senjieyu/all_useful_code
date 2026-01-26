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

mkdir -p eos_logs eos_results eos_relaxed
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
  local base tag vs oname outfile logfile

  base="$(basename "$f")"
  tag="${base%.lmp}"
  vs="$(echo "$tag" | sed -n 's/.*V\([0-9]\+\.[0-9]\+\).*/\1/p')"
  [ -z "${vs}" ] && vs="NA"

  oname="eos_relaxed/${tag}_relaxed.data"
  outfile="eos_results/${tag}.dat"
  logfile="eos_logs/${tag}.log"

  echo "[RUN] ${base}  Vscale=${vs}"

  mpirun -bootstrap lsf -np "${RANKS_PER_TASK}" "${LMP_BIN}" \
    -in "${INFILE}" \
    -log "${logfile}" \
    -var fname "${f}" \
    -var oname "${oname}" \
    -var outfile "${outfile}" \
    -var mlip_ini "${MLIP_INI}"
}

# run all points (parallel)
for f in "${FILES[@]}"; do
  wait_for_slot
  run_one "${f}" &
done
wait
echo "[INFO] All LAMMPS runs finished."

# gather results
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

