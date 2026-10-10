#!/bin/bash
#BSUB -J vasp_9
#BSUB -q 33
#BSUB -n 40
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=40]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP vasp_std"
echo "Starting relaxation (全弛豫)..."
cp 弛豫INCAR INCAR
$COMMAND_std > vasp_relax.out

echo "Starting single-point calculation (单点计算)..."
cp CONTCAR POSCAR
cp INCAR单点计算 INCAR
$COMMAND_std > vasp_scf.out

