#!/bin/bash
#BSUB -J Cu_TTT
#BSUB -q 1t88c
#BSUB -n 96
#BSUB -R "span[ptile=96]"
#BSUB -e %J.err
#BSUB -o %J.out

# --- 设置 LAMMPS 路径 (请确认路径正确) ---
LMP_EXEC="/work/phy-luoyf/apps/lammps2025/lammps-22Jul2025/build/lmp"

# 创建文件夹
mkdir -p datattt
mkdir -p relaxed_ttt

echo "==================================="
echo "Starting TTT (Full Relax) Batch..."

for vasp_file in struttt/*; do
    # 确保是文件
    [ -f "$vasp_file" ] || continue
    
    filename=$(basename "$vasp_file")
    data_file="datattt/${filename}.data"
    result_file="relaxed_ttt/${filename}_relaxed.data"

    # 1. 转换格式
    python vasp2data.py "$vasp_file" "$data_file"

    # 2. 运行 LAMMPS
    echo "Processing: $filename"
    $LMP_EXEC -in in.relax_ttt -var fname "$data_file" -var oname "$result_file" > /dev/null
done

echo "TTT Batch Finished."
echo "==================================="

