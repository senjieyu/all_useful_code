#!/bin/bash
#BSUB -J sus2_relax_sp
#BSUB -q 33
#BSUB -n 40
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=40]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
# export PATH=/work/phy-huangj/app/lammps/bin:$PATH
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP /work/phy-luoyf/apps/lammps_sus2/interface-lammps-mlip-2/lmp_intel_cpu_intelmpi"

echo "Starting LAMMPS SUS2 relaxation and single point calculation..."
$COMMAND_std -in in.sus2 > lammps.out
