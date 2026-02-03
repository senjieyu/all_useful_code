import os
import sys
import glob
import numpy as np

def get_elastic_tensor_from_outcar(outcar_path):
    """
    解析 OUTCAR，提取 TOTAL ELASTIC MODULI (kBar) 矩阵。
    返回单位为 GPa 的 6x6 numpy 数组。
    """
    if not os.path.exists(outcar_path):
        return None
        
    c_matrix = []
    found = False
    
    with open(outcar_path, 'r') as f:
        lines = f.readlines()
        
    for i, line in enumerate(lines):
        if "TOTAL ELASTIC MODULI (kBar)" in line:
            #通常下面紧接着是表头和分隔线，数据大约在发现这行后的第 3 行开始
            # VASP 格式:
            #  TOTAL ELASTIC MODULI (kBar)
            #  Direction    XX          YY    ...
            #  ----------------------------------
            #    XX      1306.40 ...
            start_idx = i + 3
            # 读取接下来的 6 行
            try:
                for j in range(6):
                    row_line = lines[start_idx + j].strip()
                    parts = row_line.split()
                    # 第一列是标签 (XX, YY 等)，后面 6 列是数据
                    # 注意：有时候 VASP 输出格式可能会连在一起，split()通常能处理空格
                    vals = [float(x) for x in parts[1:7]]
                    c_matrix.append(vals)
                found = True
                break
            except (ValueError, IndexError):
                return None

    if not found or len(c_matrix) != 6:
        return None
    
    # 转换为 numpy 数组并从 kBar 转换为 GPa (1 GPa = 10 kBar)
    c_tensor_gpa = np.array(c_matrix) / 10.0
    return c_tensor_gpa

def calculate_moduli(c_matrix):
    """
    根据 6x6 刚度矩阵 (GPa) 计算 Voigt-Reuss-Hill 模量。
    """
    # 1. Voigt bounds (Upper) using Stiffness Matrix C
    # Kv = ( (C11+C22+C33) + 2(C12+C23+C13) ) / 9
    kv = ((c_matrix[0,0] + c_matrix[1,1] + c_matrix[2,2]) + 
          2 * (c_matrix[0,1] + c_matrix[1,2] + c_matrix[0,2])) / 9.0
    
    # Gv = ( (C11+C22+C33) - (C12+C23+C13) + 3(C44+C55+C66) ) / 15
    gv = ((c_matrix[0,0] + c_matrix[1,1] + c_matrix[2,2]) - 
          (c_matrix[0,1] + c_matrix[1,2] + c_matrix[0,2]) + 
          3 * (c_matrix[3,3] + c_matrix[4,4] + c_matrix[5,5])) / 15.0
          
    # 2. Reuss bounds (Lower) using Compliance Matrix S = C^-1
    try:
        s_matrix = np.linalg.inv(c_matrix)
    except np.linalg.LinAlgError:
        return None # 矩阵奇异，无法求逆
        
    # Kr = 1 / [ (S11+S22+S33) + 2(S12+S23+S13) ]
    s_sum_k = (s_matrix[0,0] + s_matrix[1,1] + s_matrix[2,2]) + \
              2 * (s_matrix[0,1] + s_matrix[1,2] + s_matrix[0,2])
    kr = 1.0 / s_sum_k
    
    # Gr = 15 / [ 4(S11+S22+S33) - 4(S12+S23+S13) + 3(S44+S55+S66) ]
    s_sum_g = 4 * (s_matrix[0,0] + s_matrix[1,1] + s_matrix[2,2]) - \
              4 * (s_matrix[0,1] + s_matrix[1,2] + s_matrix[0,2]) + \
              3 * (s_matrix[3,3] + s_matrix[4,4] + s_matrix[5,5])
    gr = 15.0 / s_sum_g
    
    # 3. Hill average
    k_vrh = (kv + kr) / 2.0
    g_vrh = (gv + gr) / 2.0
    
    # 4. Poisson's Ratio (using VRH values)
    # nu = (3K - 2G) / (2(3K + G))
    if (3*k_vrh + g_vrh) != 0:
        nu = (3*k_vrh - 2*g_vrh) / (2 * (3*k_vrh + g_vrh))
    else:
        nu = 0.0
        
    # 5. Universal Anisotropy Index
    # Au = 5*(Gv/Gr) + (Kv/Kr) - 6
    if gr != 0 and kr != 0:
        anisotropy = 5 * (gv / gr) + (kv / kr) - 6.0
    else:
        anisotropy = 0.0
        
    return {
        "Kv": kv, "Kr": kr, "Kvrh": k_vrh,
        "Gv": gv, "Gr": gr, "Gvrh": g_vrh,
        "nu": nu, "anisotropy": anisotropy
    }

def format_output(c_matrix, props):
    """格式化输出字符串"""
    lines = []
    lines.append("Stiffness Tensor Cij (GPa)")
    
    # Format Matrix
    for row in c_matrix:
        # 使用 10.2f 格式对齐
        line_str = "".join([f"{val:10.2f} " for val in row])
        lines.append(line_str)
        
    if props:
        lines.append(f"Bulk Modulus Kv:   {props['Kv']:.2f} GPa")
        lines.append(f"Bulk Modulus Kr:   {props['Kr']:.2f} GPa")
        lines.append(f"Bulk Modulus Kvrh: {props['Kvrh']:.2f} GPa")
        lines.append(f"Shear Modulus Gv:   {props['Gv']:.2f} GPa")
        lines.append(f"Shear Modulus Gr:   {props['Gr']:.2f} GPa")
        lines.append(f"Shear Modulus Gvrh: {props['Gvrh']:.2f} GPa")
        lines.append(f"Elastic Anisotropy: {props['anisotropy']:.2f}")
        lines.append(f"Poisson's Ratio:    {props['nu']:.2f}")
    else:
        lines.append("Calculation failed (Singular Matrix?)")
        
    return "\n".join(lines)

def main():
    if len(sys.argv) < 3:
        print("Usage: python collect_vasp_elastic.py <work_dir> <final_out_dir>")
        sys.exit(1)
        
    work_dir = sys.argv[1]
    final_out_dir = sys.argv[2]
    
    if not os.path.exists(work_dir):
        print(f"Work directory {work_dir} not found.")
        sys.exit(1)
        
    os.makedirs(final_out_dir, exist_ok=True)
    
    struct_dirs = sorted(glob.glob(os.path.join(work_dir, "struct_*")))
    
    count = 0
    for s_dir in struct_dirs:
        s_name = os.path.basename(s_dir) # struct_0
        outcar = os.path.join(s_dir, "OUTCAR")
        
        c_matrix_gpa = get_elastic_tensor_from_outcar(outcar)
        
        if c_matrix_gpa is not None:
            props = calculate_moduli(c_matrix_gpa)
            result_text = format_output(c_matrix_gpa, props)
            
            out_file = os.path.join(final_out_dir, f"{s_name}_elastic.txt")
            with open(out_file, 'w') as f:
                f.write(f"Results for {s_name}\n")
                f.write("-" * 30 + "\n")
                f.write(result_text)
            count += 1
        else:
            print(f"Warning: OUTCAR not found or elastic data missing in {s_name}")
            
    print(f"Collected and processed results for {count} structures to {final_out_dir}")

if __name__ == "__main__":
    main()

