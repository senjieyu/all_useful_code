#!/bin/bash
#BSUB -J GaAs_VASP
#BSUB -q 1t81c
#BSUB -n 16
#BSUB -W 1:00
#BSUB -o VASP_%J.out
#BSUB -e VASP_%J.err
#BSUB -R "span[ptile=16]"

# ===== 基本环境设置 =====
export VASP_ROOT=/share/apps/vasp
export VASP_BIN=$VASP_ROOT/vasp.6.3.0_IceLake/bin
export PATH=$VASP_BIN:$PATH

# oneAPI 库路径（根据 ldd 结果）
export LD_LIBRARY_PATH=/share/intel/oneapi/mkl/2022.2.0/lib/intel64:/share/intel/oneapi/mpi/2021.6.0/lib/release:/share/intel/oneapi/mpi/2021.6.0/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=1
ulimit -s unlimited

# ===== 运行 VASP =====
VASP_EXE=$VASP_BIN/vasp_std
NP=${LSB_DJOB_NUMPROC:-16}

echo "Running on $NP cores..."
mpirun -bootstrap lsf -envall -np $NP $VASP_EXE

# ===== 结束提示 =====
if [ $? -eq 0 ]; then
  echo "✅  VASP completed successfully at $(date)"
else
  echo "❌  VASP failed at $(date)"
fi

