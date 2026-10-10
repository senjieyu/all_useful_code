#!/bin/bash
#BSUB -J Al2_d02
#BSUB -q 33
#BSUB -n 40
#BSUB -R "span[ptile=40]"
#BSUB -o %J.out
#BSUB -e %J.err

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH
mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP vasp_std