import sys
import os
import numpy as np

# --- 铜原子质量 ---
ATOM_MASS = 63.546

def convert_poscar_to_data(poscar_file, output_file):
    if not os.path.exists(poscar_file):
        print(f"Error: {poscar_file} not found.")
        return

    try:
        with open(poscar_file, 'r') as f:
            lines = f.readlines()
            
        if len(lines) < 7:
            print(f"Warning: {os.path.basename(poscar_file)} looks too short, skipping.")
            return

        # 1. 读取缩放系数
        try:
            scale = float(lines[1].strip())
        except ValueError:
            print(f"Error reading scale in {poscar_file}")
            return

        # 2. 读取晶格矩阵
        lattice = []
        for i in range(2, 5):
            lattice.append([float(x) for x in lines[i].split()])
        lattice = np.array(lattice) * scale

        # --- 核心修正：计算 LAMMPS 斜盒子参数 (Triclinic Box) ---
        a = lattice[0]
        b = lattice[1]
        c = lattice[2]
        
        lx = np.linalg.norm(a)
        unit_a = a / lx if lx > 1e-8 else np.array([1.0, 0.0, 0.0])
        
        xy = np.dot(b, unit_a)
        ly = np.sqrt(max(0.0, np.linalg.norm(b)**2 - xy**2))
        
        xz = np.dot(c, unit_a)
        yz = (np.dot(b, c) - xy*xz) / ly if ly > 1e-8 else 0.0
        lz = np.sqrt(max(0.0, np.linalg.norm(c)**2 - xz**2 - yz**2))

        # 3. 确定原子数行 (兼容 VASP 4/5)
        try:
            [int(x) for x in lines[5].split()]
            natoms_line_idx = 5
            start_idx = 6
        except ValueError:
            natoms_line_idx = 6
            start_idx = 7
        
        natoms_list = [int(x) for x in lines[natoms_line_idx].split()]
        total_atoms = sum(natoms_list)
        
        # 4. 检测 Selective Dynamics
        selective = False
        if lines[start_idx].strip().upper().startswith('S'):
            selective = True
            start_idx += 1
        
        coord_type = lines[start_idx].strip().lower()
        start_idx += 1

        atoms_data = []
        atom_id = 1
        current_line = start_idx
        
        # 5. 读取坐标
        for n in natoms_list:
            for _ in range(n):
                if current_line >= len(lines): break
                line_content = lines[current_line].split()
                
                # --- 类型分配逻辑 ---
                atom_type = 1 
                if selective and len(line_content) >= 6:
                    flags = [f.upper() for f in line_content[3:6]]
                    if flags == ['F', 'F', 'F']:
                        atom_type = 2 # FFF
                    elif flags == ['F', 'F', 'T']:
                        atom_type = 3 # FFT
                
                # 坐标转换：Direct -> Cartesian
                vals = [float(v) for v in line_content[0:3]]
                if coord_type.startswith('d'): # Direct
                    vec = np.array(vals)
                    cart_coords = np.dot(vec, lattice)
                    x, y, z = cart_coords
                else: # Cartesian
                    x = vals[0] * scale
                    y = vals[1] * scale
                    z = vals[2] * scale

                atoms_data.append((atom_id, atom_type, x, y, z))
                atom_id += 1
                current_line += 1

        # 6. 写入 LAMMPS Data
        with open(output_file, 'w') as f:
            f.write(f"Converted from {os.path.basename(poscar_file)} (Triclinic/Cu)\n\n")
            f.write(f"{total_atoms} atoms\n")
            f.write("3 atom types\n")
            
            f.write(f"0.0 {lx:.6f} xlo xhi\n")
            f.write(f"0.0 {ly:.6f} ylo yhi\n")
            f.write(f"0.0 {lz:.6f} zlo zhi\n")
            f.write(f"{xy:.6f} {xz:.6f} {yz:.6f} xy xz yz\n")
            
            f.write("\nMasses\n\n")
            f.write(f"1 {ATOM_MASS}\n")
            f.write(f"2 {ATOM_MASS}\n")
            f.write(f"3 {ATOM_MASS}\n")
            
            f.write("\nAtoms # atomic\n\n")
            for data in atoms_data:
                f.write(f"{data[0]} {data[1]} {data[2]:.6f} {data[3]:.6f} {data[4]:.6f}\n")
        
        print(f"Converted: {os.path.basename(poscar_file)} -> {os.path.basename(output_file)}")

    except Exception as e:
        print(f"Error converting {os.path.basename(poscar_file)}: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python vasp2data.py <input_folder_or_file> <output_folder_or_file>")
        print("Example: python vasp2data.py struttt datattt")
        sys.exit(1)

    input_arg = sys.argv[1]
    output_arg = sys.argv[2]

    # --- 文件夹模式 (Batch Processing) ---
    if os.path.isdir(input_arg):
        # 如果输出目录不存在，自动创建
        if not os.path.exists(output_arg):
            os.makedirs(output_arg)
            print(f"Created output directory: {output_arg}")
        
        print(f"Processing directory: {input_arg} -> {output_arg}")
        
        file_list = sorted(os.listdir(input_arg))
        count = 0
        for filename in file_list:
            src_file = os.path.join(input_arg, filename)
            
            # 跳过子文件夹和隐藏文件
            if not os.path.isfile(src_file) or filename.startswith('.'):
                continue
            
            # 构造输出文件名: 去掉后缀，加上 .data
            # 例如: structure.vasp -> structure.data
            base_name = os.path.splitext(filename)[0]
            dst_file = os.path.join(output_arg, base_name + ".data")
            
            convert_poscar_to_data(src_file, dst_file)
            count += 1
            
        print(f"Done! Processed {count} files.")

    # --- 单文件模式 (Single File) ---
    else:
        # 如果输出路径包含文件夹且不存在，尝试创建
        out_dir = os.path.dirname(output_arg)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)
        convert_poscar_to_data(input_arg, output_arg)

