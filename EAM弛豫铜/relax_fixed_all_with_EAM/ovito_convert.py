import sys
import os
import glob
from ovito.io import import_file, export_file
from ovito.data import ParticleProperty

# 检查是否提供了输入文件夹
if len(sys.argv) != 3:
    print("Usage: python ovito_convert.py <input_folder> <output_folder>")
    print("Example: python ovito_convert.py strufft datafft")
    sys.exit(1)

input_folder = sys.argv[1]
output_folder = sys.argv[2]

# 确保输出目录存在
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# 获取所有 .vasp 文件
files = glob.glob(os.path.join(input_folder, "*.vasp"))
if not files:
    print(f"No .vasp files found in {input_folder}")
    sys.exit(1)

print(f"Found {len(files)} files in {input_folder}...")

for file_path in files:
    filename = os.path.basename(file_path)
    output_path = os.path.join(output_folder, filename.replace('.vasp', '.data'))
    
    print(f"Converting: {filename} -> {os.path.basename(output_path)}")
    
    # 1. 导入 VASP 文件
    pipeline = import_file(file_path)
    
    # 2. 定义修改器函数：将 Selective Dynamics 映射为 Particle Type
    def modify_particle_types(frame, data):
        # 检查是否存在 Selective Dynamics 属性
        # VASP 导入后，Selective Dynamics 通常存储为 "Selective Dynamics" 属性 (0=False, 1=True ? 或相反，需验证)
        # OVITO 中: Selective Dynamics 属性通常是 (Int, Int, Int)，1代表T(移动)，0代表F(固定)
        
        if 'Selective Dynamics' in data.particles:
            sel_dyn = data.particles['Selective Dynamics']
            types = data.particles_.particle_types_ # 获取类型数组 (可写)
            
            # 这里的逻辑需要根据 OVITO 版本确认，通常 1=T (Active), 0=F (Fixed)
            # 我们的目标:
            # T T T (1 1 1) -> Type 1
            # F F F (0 0 0) -> Type 2
            # F F T (0 0 1) -> Type 3
            
            # 遍历所有粒子 (使用 numpy 加速更佳，但为了清晰用循环)
            # 创建新的类型数组
            new_types = data.particles.particle_types[:] 
            
            for i in range(data.particles.count):
                flags = sel_dyn[i] # [x, y, z]
                
                if flags[0] == 0 and flags[1] == 0 and flags[2] == 0:
                    new_types[i] = 2 # FFF
                elif flags[0] == 0 and flags[1] == 0 and flags[2] == 1:
                    new_types[i] = 3 # FFT
                else:
                    new_types[i] = 1 # TTT 或其他
            
            # 更新类型
            data.particles_.particle_types_ = new_types
        else:
            # 如果没有 Selective Dynamics 属性，全设为 Type 1
            pass

    # 3. 应用修改器
    pipeline.modifiers.append(modify_particle_types)
    
    # 4. 计算并导出
    data = pipeline.compute()
    
    # 导出为 LAMMPS Data 格式
    # atom_style atomic 对应 id type x y z
    export_file(
        data, 
        output_path, 
        "lammps/data", 
        atom_style="atomic"
    )

print("Batch conversion finished.")

