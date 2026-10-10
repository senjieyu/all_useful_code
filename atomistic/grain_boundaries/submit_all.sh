#!/bin/bash
# 一键准备和提交三种方法的作业 (请在集群上执行)

echo ">>> Submitting VASP job..."
cd VASP
bsub < sub.sh
cd ..

echo ">>> Converting POSCAR to LAMMPS data format..."
python3 vasp2lammps.py POSCAR LAMMPS_EAM/data.lammps
python3 vasp2lammps.py POSCAR LAMMPS_SUS2/data.lammps

echo ">>> Submitting LAMMPS EAM job..."
cd LAMMPS_EAM
bsub < sub.sh
cd ..

echo ">>> Submitting LAMMPS SUS2 job..."
cd LAMMPS_SUS2
bsub < sub.sh
cd ..

echo "All jobs submitted!"
