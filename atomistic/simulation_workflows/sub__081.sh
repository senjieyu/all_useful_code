#!/bin/bash
#BSUB -J pred_153
#BSUB -q 33
#BSUB -n 40
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=40]"

# Set up environment
# Note: Ensure python, ase, numpy, matplotlib, sklearn are available
# module load python/3.9 # Example: Uncomment if needed
# source /path/to/your/venv/bin/activate # Example: Uncomment if needed

hostfile=$LSB_DJOB_HOSTFILE
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD

# Define LAMMPS executable with MLIP support
# Make sure this path is correct on your server
LMP_EXEC="/work/phy-huangj/app/interface-lammps-mlip-v1.2-2/lmp_intel_cpu_intelmpi"

# Construct the full command for LAMMPS
# This variable will be used by the python script
export LMP_BIN="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP $LMP_EXEC"

echo "Running calculation with command: $LMP_BIN"

# Run the python script
python xyz2lammps_and_eval.py
