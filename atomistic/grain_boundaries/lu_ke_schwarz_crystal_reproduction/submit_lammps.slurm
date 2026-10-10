#!/bin/bash
#SBATCH --job-name=test                 # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=your_email@example.com  # Optional notification address
#SBATCH --nodes=1                     # Maximum number of nodes to be allocated
#SBATCH --ntasks=12			# Run on 96 cores
#SBATCH --partition=P                   # Partitions name
#SBATCH --output=test_%j.log            # Standard output and error log

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

# Adapt module names to your cluster.
# source /path/to/oneapi/setvars.sh
# module load lammps
# module load cuda
export UCX_TLS=ib
#--------------------------Run on CPU only-------------------------------
#mpirun -np $SLURM_NTASKS lmp_gpu -in latticeVCoNi.in

#---------------------Accelerated by GPU-----------------------------------
mpirun -np $SLURM_NTASKS lmp_gpu -sf gpu -pk gpu 1 -in melt.in

#---------------------Accelerated by GPU with neigh no---------------------
#mpirun -np $SLURM_NTASKS lmp_gpu -sf gpu -pk gpu 1 neigh no -in in.lammps.in
