#!/bin/bash
#BSUB -J Cu_POSCAR_Cu111_twin_singleplane_d1_7_d2_12_d3_5_1100.vasp_relaxed_relax
#BSUB -q v3-64
#BSUB -n 26
#BSUB -R "span[ptile=24]"
#BSUB -e relax_%J.err
#BSUB -o relax_%J.out

set -euo pipefail

hostfile="$LSB_DJOB_HOSTFILE"
if [ -z "$hostfile" ]; then
    NP=4
    CMD="mpirun -np $NP"
else
    NP=$(wc -l < "$hostfile")
    CMD="mpirun -machinefile $hostfile -np $NP"
fi

cd "$LS_SUBCWD"
module load mpi/2021.6.0 compiler/2022.1.0 mkl/2022.2.0

echo "=========================================="
echo "Relax Job Start: $(date)"
echo "Structure: POSCAR_Cu111_twin_singleplane_d1_7_d2_12_d3_5_1100.vasp_relaxed"
echo "=========================================="

$CMD "/work/phy-luoyf/apps/lammps2025/lammps-22Jul2025/build/lmp" -in "relax.in" -var struct_file "POSCAR_Cu111_twin_singleplane_d1_7_d2_12_d3_5_1100.vasp_relaxed.data" > log.relax

if [ ! -f relaxed.data ]; then
    echo "Error: relaxed.data not found after relaxation."
    exit 1
fi

cp relaxed.data ../stru_relaxed/POSCAR_Cu111_twin_singleplane_d1_7_d2_12_d3_5_1100.vasp_relaxed_relaxed.data

echo "=========================================="
echo "Relax Job Done: $(date)"
echo "=========================================="
