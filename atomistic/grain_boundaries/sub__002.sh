#!/bin/bash
#BSUB -J eam_1
#BSUB -q 33
#BSUB -n 40
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=40]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
# 若有需要，请取消注释并修改为你的LAMMPS路径
# export PATH=/work/phy-huangj/app/lammps/bin:$PATH
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP lmp"

echo "Starting LAMMPS EAM calculation..."
$COMMAND_std -in in.eam > lammps.out
