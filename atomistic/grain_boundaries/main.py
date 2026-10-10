# Main Entry Point
# 程序入口文件

import numpy as np
import setting
from gb_code import csl_generator as csl
from gb_code import gb_generator as gb

def step1_search_sigma():
    """
    Step 1: 寻找可能的 Sigma 值
    根据 setting.py 中的 AXIS (旋转轴) 和 LIMIT_SIGMA (限制) 进行搜索。
    """
    print(f"=== Step 1: Searching CSL for Axis: {setting.AXIS} ===")
    uvw = np.array(setting.AXIS)
    csl.print_list(uvw, setting.LIMIT_SIGMA)
    print("\n[提示] 请从上方列表中选择一个 Sigma 值，并更新 setting.py 中的 'SIGMA' 参数。")

def step2_search_planes():
    """
    Step 2: 寻找可能的晶界平面
    根据 setting.py 中的 SIGMA, BASIS 进行搜索。
    """
    print(f"=== Step 2: Searching GB Planes for Sigma: {setting.SIGMA} ===")
    uvw = np.array(setting.AXIS)
    sigma = setting.SIGMA
    
    # 获取 theta, m, n
    thetas = csl.get_theta_m_n_list(uvw, sigma)
    if not thetas:
        print(f"Error: No valid CSL found for Sigma {sigma}")
        return
    
    theta, m, n = thetas[0]
    
    print("-----------------------------------------------------------")
    print(" GB1                   GB2                   Type          Number of Atoms")
    csl.print_list_GB_Planes(uvw, setting.BASIS, m, n, setting.LIMIT_GB_PLANE)
    print("\n[提示] 请从 'GB1' 列中选择一个晶界平面，并更新 setting.py 中的 'GB_PLANE' 参数。")

def step3_generate_gb():
    """
    Step 3: 生成晶界原子结构
    根据 setting.py 中的所有参数生成最终的结构文件。
    """
    print(f"=== Step 3: Generating GB Structure for Plane: {setting.GB_PLANE} ===")
    
    # 准备参数
    uvw = np.array(setting.AXIS)
    sigma = setting.SIGMA
    thetas = csl.get_theta_m_n_list(uvw, sigma)
    if not thetas:
        print(f"Error: No valid CSL found for Sigma {sigma}")
        return
    theta, m, n = thetas[0]
    
    # 初始化 GB 对象
    my_gb = gb.GB_character()
    
    # 解析参数
    print("Parsing parameters...")
    my_gb.ParseGB(
        axis=uvw,
        basis=setting.BASIS,
        LatP=setting.LATTICE_PARAMETER,
        m=m,
        n=n,
        gb=setting.GB_PLANE,
        atom_type=setting.ATOM_TYPE,
        auto_lat=setting.AUTO_LATTICE_PARAMETER
    )
    
    # 生成双晶结构
    print("Generating Bicrystal atoms...")
    my_gb.CSL_Bicrystal_Atom_generator()
    
    # 写入文件
    print(f"Writing output to {setting.FILE_TYPE} format...")
    my_gb.WriteGB(
        overlap=setting.OVERLAP_DISTANCE,
        whichG=setting.WHICH_G,
        rigid=setting.RIGID_TRANS,
        a=setting.RIGID_A,
        b=setting.RIGID_B,
        dim1=setting.DIMENSIONS[0],
        dim2=setting.DIMENSIONS[1],
        dim3=setting.DIMENSIONS[2],
        file=setting.FILE_TYPE
    )
    print("\n[完成] 结构生成完毕。请检查生成的输出文件。")

if __name__ == "__main__":
    print("###############################################")
    print("#           GB Code Automation Tool           #")
    print("###############################################")
    print("当前配置读取自 setting.py")
    print("请取消注释下方对应的函数调用来执行特定步骤：\n")
    
    # 1. 第一步：寻找 Sigma
    # step1_search_sigma()
    
    # 2. 第二步：寻找晶界平面
    # step2_search_planes()
    
    # 3. 第三步：生成原子结构 (默认开启)
    step3_generate_gb()
