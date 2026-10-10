#!/bin/bash
#BSUB -J Cu_100_50
#BSUB -q 33
#BSUB -n 40
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=40]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH

# Link the special binary
if [ -f ../../../vasp_std_str_c ]; then
    ln -sf ../../../vasp_std_str_c .
fi

COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP ./vasp_std_str_c"

echo "Starting VASP relaxation (ISIF=3, vasp_std_str_c)..."
$COMMAND_std > vasp_relax.log
