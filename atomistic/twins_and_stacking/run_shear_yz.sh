#!/bin/bash
#BSUB -J Cu_POSCAR_Cu111_twin_singleplane_d1_7_d2_12_d3_5_1100.vasp_relaxed_yz
#BSUB -q v3-64
#BSUB -n 26
#BSUB -R "span[ptile=24]"
#BSUB -e shear_yz_%J.err
#BSUB -o shear_yz_%J.out

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
echo "Shear Job Start: $(date)"
echo "Structure: POSCAR_Cu111_twin_singleplane_d1_7_d2_12_d3_5_1100.vasp_relaxed | Mode: yz"
echo "=========================================="

mkdir -p yz
cd yz

ln -sf ../relaxed.data .
ln -sf ../Cu_mishin1.eam.alloy .

$CMD "/work/phy-luoyf/apps/lammps2025/lammps-22Jul2025/build/lmp" -in ../"tensile.in" -var struct_file relaxed.data -var shear_mode yz > log.shear_yz

python ../"analyze.py" stress_strain_yz.dat stress_strain_yz.png yield_info.txt yz

GLOBAL_RES_DIR="../../results_YZ"
cp stress_strain_yz.png "$GLOBAL_RES_DIR/POSCAR_Cu111_twin_singleplane_d1_7_d2_12_d3_5_1100.vasp_relaxed_stress_strain.png"
cp stress_strain_yz.dat "$GLOBAL_RES_DIR/POSCAR_Cu111_twin_singleplane_d1_7_d2_12_d3_5_1100.vasp_relaxed_stress_strain.dat"
cp final_structure_yz.data "$GLOBAL_RES_DIR/POSCAR_Cu111_twin_singleplane_d1_7_d2_12_d3_5_1100.vasp_relaxed_final.data"

echo "Structure: POSCAR_Cu111_twin_singleplane_d1_7_d2_12_d3_5_1100.vasp_relaxed | Mode: yz" >> "/work/phy-luoyf/WORK/cu_sample_construct_after_wyb_meeting25.12.25/26.3.29_CU_complex_example/xxx_shear/shear_summary.txt"
cat yield_info.txt >> "/work/phy-luoyf/WORK/cu_sample_construct_after_wyb_meeting25.12.25/26.3.29_CU_complex_example/xxx_shear/shear_summary.txt"
echo "" >> "/work/phy-luoyf/WORK/cu_sample_construct_after_wyb_meeting25.12.25/26.3.29_CU_complex_example/xxx_shear/shear_summary.txt"

echo "=========================================="
echo "Shear Job Done: $(date)"
echo "=========================================="
