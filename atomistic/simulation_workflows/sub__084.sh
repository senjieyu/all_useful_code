#!/bin/bash
#BSUB -J scf_job
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

echo "Starting single-point calculation (单点计算)..."
$COMMAND_std > vasp_scf.out

echo "Calculating average energy..."
python3 calc_energy.py > energy_result.txt
cat energy_result.txt
