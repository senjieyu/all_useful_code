#!/bin/bash
#BSUB -J Cu_Relax_Batch
#BSUB -q 1t88c
#BSUB -n 96
#BSUB -R "span[ptile=96]"
#BSUB -e %J.err
#BSUB -o %J.out

# --- 设置路径 ---
# LAMMPS 可执行文件绝对路径
LMP_EXEC="/work/phy-luoyf/apps/lammps2025/lammps-22Jul2025/build/lmp"
# Python 脚本路径 (假设就在当前目录)
PY_SCRIPT="vasp2data.py"

# --- 准备文件夹 ---
# 创建存放转换后 data 文件的目录
if [ ! -d "stru_data" ]; then
    mkdir stru_data
fi

# 创建存放优化结果的目录
if [ ! -d "relaxed" ]; then
    mkdir relaxed
fi

echo "========================================"
echo "Job started at $(date)"
echo "Step 1: Converting VASP files to LAMMPS Data format..."

# --- 循环处理 stru/ 下的所有文件 ---
for vasp_file in stru/*; do
    
    # 跳过非文件项
    [ -f "$vasp_file" ] || continue
    
    # 获取文件名 (例如 1D_main_000.vasp)
    filename=$(basename "$vasp_file")
    # 构造 data 文件名 (例如 1D_main_000.data)
    data_file="stru_data/${filename%.*}.data"
    # 构造输出结果文件名
    relaxed_file="relaxed/${filename%.*}_relaxed.data"

    # 1. 执行 Python 转换
    # 自动识别 TTT/FFF/FFT 并生成含 Type 1/2/3 的 data 文件
    python $PY_SCRIPT "$vasp_file" "$data_file"
    
    if [ ! -f "$data_file" ]; then
        echo "Warning: Failed to convert $filename"
        continue
    fi

    echo "----------------------------------------"
    echo "Processing: $filename"
    echo "  Mode: Relaxing (Type 2 fixed, Type 3 z-only)"

    # 2. 运行 LAMMPS
    # 注意：小体系强烈建议单核运行，不要用 mpirun
    $LMP_EXEC -in in.relax \
              -var fname "$data_file" \
              -var oname "$relaxed_file" \
              > /dev/null

    if [ $? -eq 0 ]; then
        echo "  Done. Result: $relaxed_file"
    else
        echo "  Error: LAMMPS failed for $filename"
    fi

done

echo "========================================"
echo "All tasks finished at $(date)"

