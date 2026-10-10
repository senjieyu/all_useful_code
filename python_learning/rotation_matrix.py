import numpy as np


def normalize(v):
    """归一化向量"""
    norm = np.linalg.norm(v)
    if norm == 0:
        raise ValueError("向量模长不能为0")
    return np.array(v) / norm


def get_basis_matrix(primary_axis_name, primary_vec, secondary_vec):
    """
    构建基变换矩阵
    :param primary_axis_name: 主轴名称 ('x', 'y', 'z')
    :param primary_vec: 主轴对应的晶向向量 (如 [1, 1, 1])
    :param secondary_vec: 辅助向量，用于确定平面 (如 [-1, 1, 0])
    :return: 3x3 旋转矩阵, (x_vec, y_vec, z_vec) 计算出的三个正交向量
    """
    # 1. 归一化主轴
    u_primary = normalize(primary_vec)

    # 2. 处理辅助轴 (Gram-Schmidt 正交化)
    # 我们需要一个垂直于主轴的向量。
    # 先归一化辅助向量
    v_temp = normalize(secondary_vec)

    # 计算辅助轴在主轴上的投影，并减去它，得到垂直分量
    # v_orth = v - (v . u) * u
    proj = np.dot(v_temp, u_primary)
    u_secondary = v_temp - proj * u_primary

    # 如果辅助向量和主轴平行，通过叉乘自动找一个
    if np.linalg.norm(u_secondary) < 1e-6:
        print("警告: 辅助向量与主轴平行，将自动选择参考向量。")
        # 尝试用 [1,0,0] 或 [0,1,0]
        temp = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(temp, u_primary)) > 0.9:
            temp = np.array([0.0, 1.0, 0.0])
        u_secondary = temp - np.dot(temp, u_primary) * u_primary

    u_secondary = normalize(u_secondary)

    # 3. 计算第三个轴 (叉乘)
    # 右手定则: X x Y = Z, Y x Z = X, Z x X = Y
    u_tertiary = np.cross(u_primary, u_secondary)  # 默认顺序，后面调整

    # 4. 根据用户指定的主轴是 X, Y 还是 Z 来分配向量
    name = primary_axis_name.lower()

    if name == 'x':
        # 主轴是 X，辅助轴通常定义 Y (或者 Z)
        # 这里假设辅助向量意图指向 Y 方向
        new_x = u_primary
        new_y = u_secondary  # 辅助轴修正后作为 Y
        new_z = np.cross(new_x, new_y)
    elif name == 'y':
        # 主轴是 Y (用户示例情况)，辅助轴假设指向 Z
        new_y = u_primary
        new_z = u_secondary  # 辅助轴修正后作为 Z
        new_x = np.cross(new_y, new_z)  # Y x Z = X
    elif name == 'z':
        # 主轴是 Z，辅助轴假设指向 X
        new_z = u_primary
        new_x = u_secondary
        new_y = np.cross(new_z, new_x)
    else:
        raise ValueError("主轴必须是 x, y 或 z")

    # 5. 构建矩阵 (以列向量形式堆叠)
    # 矩阵 R 的列就是新坐标系的基向量在旧坐标系中的表示
    R = np.column_stack((new_x, new_y, new_z))

    return R, (new_x, new_y, new_z)


def main():
    print("--- 晶体取向矩阵计算器 (Basis Alignment) ---")
    print("目标: 将 XYZ 坐标轴对准特定的晶向 (Miller Indices)")

    try:
        # 1. 获取主轴输入 (用户想让 Y 轴变成 [1 1 1])
        axis_name = input("请输入主轴名称 (x/y/z) [默认 y]: ").strip().lower() or 'y'

        vec_str = input(f"请输入主轴 {axis_name.upper()} 的方向向量 (例如 '1 1 1'): ")
        primary_vec = np.array([float(x) for x in vec_str.split()])

        # 2. 获取辅助轴输入 (为了确定唯一旋转，必须再给一个向量)
        print("\n需要辅助向量来固定旋转角度 (确定另外两个轴的方向)。")
        print(f"如果你想让结果出现 [-1 1 0]，请在这里输入它。")
        aux_str = input("请输入辅助方向向量 (例如 '-1 1 0'): ")
        secondary_vec = np.array([float(x) for x in aux_str.split()])

        # 3. 计算
        R, vectors = get_basis_matrix(axis_name, primary_vec, secondary_vec)

        # 4. 输出结果
        print("\n" + "=" * 30)
        print(f"计算出的旋转矩阵 R (需乘的矩阵):")
        with np.printoptions(precision=4, suppress=True):
            print(R)

        print("\n新坐标轴对应的单位向量 (列向量):")
        print(f"New X: {vectors[0]}")
        print(f"New Y: {vectors[1]}")
        print(f"New Z: {vectors[2]}")

        print("-" * 30)
        # 尝试还原为整数比例 (帮助用户验证是否是 [1 1 -2] 这种形式)
        print("归一化向量的整数比例猜测 (仅供参考):")
        for label, v in zip(['X', 'Y', 'Z'], vectors):
            # 简单算法：除以最小的非零绝对值，看是否接近整数
            min_val = np.min(np.abs(v[np.abs(v) > 1e-4]))
            scaled = v / min_val
            print(f"New {label} 方向 ≈ {np.round(scaled, 1)}")

    except Exception as e:
        print(f"错误: {e}")


if __name__ == "__main__":
    main()
