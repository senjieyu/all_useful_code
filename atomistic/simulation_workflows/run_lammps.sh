#!/bin/bash
#BSUB -J Cu_bulk_uniaxial
#BSUB -q 1t88c
#BSUB -n 96
#BSUB -R "span[ptile=96]"
#BSUB -e %J.err
#BSUB -o %J.out

# ---- environment ----
hostfile="$LSB_DJOB_HOSTFILE"
NP=$(wc -l < "$hostfile")
cd "$LS_SUBCWD"

module load mpi/2021.6.0 compiler/2022.1.0 mkl/2022.2.0

# === your LAMMPS (absolute path, no ~)
LAMMPS_EXE="/work/phy-luoyf/apps/lammps2025/lammps-22Jul2025/build/lmp"

# Input & potential
INPUT="in_cu_nanocube_bulk.in"
POT="Cu_mishin1.eam.alloy"

OUTDIR="./output_bulk"
mkdir -p "$OUTDIR"

export OMP_NUM_THREADS=1

# ---- sanity checks ----
[[ -x "$LAMMPS_EXE" ]] || { echo "LAMMPS binary NOT found"; exit 1; }
[[ -f "$INPUT" ]] || { echo "Input file NOT found"; exit 2; }
grep -q "$POT" "$INPUT" && [[ ! -f "$POT" ]] && { echo "Potential missing"; exit 3; }

echo "LAMMPS: $LAMMPS_EXE"
echo "NP=$NP"
echo "Input: $INPUT"
echo "OutputDir: $OUTDIR"

# ---- run ----
start_time=$(date +%s.%N)

mpirun -machinefile "$hostfile" -np "$NP" "$LAMMPS_EXE" -in "$INPUT" \
  | tee "$OUTDIR/lammps.log"

rc=${PIPESTATUS[0]}

end_time=$(date +%s.%N)
runtime=$(echo "$end_time - $start_time" | bc)
echo "total_runtime: $runtime s" | tee -a "$OUTDIR/time.txt"

[[ $rc -ne 0 ]] && { echo "LAMMPS exited with non-zero status: $rc"; exit $rc; }

echo "LAMMPS job finished OK."

