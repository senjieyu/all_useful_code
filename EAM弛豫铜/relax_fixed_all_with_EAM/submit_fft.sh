#!/bin/bash
#BSUB -J Cu_FFT
#BSUB -q 1t88c
#BSUB -n 96
#BSUB -R "span[ptile=96]"
#BSUB -e %J.err
#BSUB -o %J.out

# --- 设置 LAMMPS 路径 ---
LMP_EXEC="/work/phy-luoyf/apps/lammps2025/lammps-22Jul2025/build/lmp"

# 创建文件夹
mkdir -p datafft
mkdir -p relaxed_fft

echo "==================================="
echo "Starting FFT (Constrained Relax) Batch..."

for vasp_file in strufft/*; do
    # 确保是文件
    [ -f "$vasp_file" ] || continue
    
    filename=$(basename "$vasp_file")
    data_file="datafft/${filename}.data"
    result_file="relaxed_fft/${filename}_relaxed.data"

    # 1. 转换格式 (Python脚本会自动识别F/T标记并设为Type 2或3)
    python vasp2data.py "$vasp_file" "$data_file"

    # 2. 运行 LAMMPS
    echo "Processing: $filename"
    $LMP_EXEC -in in.relax_fft -var fname "$data_file" -var oname "$result_file" > /dev/null
done

echo "FFT Batch Finished."
echo "==================================="

