#!/bin/bash
#BSUB -J relax_job
#BSUB -q 1t88c
#BSUB -n 96
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=96]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP vasp_std"

echo "Starting relaxation (全弛豫)..."
$COMMAND_std > vasp_relax.out
echo "Relaxation finished. Please use the CONTCAR file for subsequent calculations."
