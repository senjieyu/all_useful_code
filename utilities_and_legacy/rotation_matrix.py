import numpy as np


def get_rotation_matrix(axis, angle_degrees):
    """
    计算绕任意轴旋转的旋转矩阵 (罗德里格旋转公式)

    参数:
    axis (list or np.array): 旋转轴向量 [x, y, z]
    angle_degrees (float): 旋转角度 (度)

    返回:
    np.array: 3x3 旋转矩阵
    """
    # 1. 将轴向量转换为 numpy 数组
    axis = np.array(axis, dtype=float)

    # 2. 归一化旋转轴 (确保它是单位向量)
    norm = np.linalg.norm(axis)
    if norm == 0:
        raise ValueError("旋转轴不能为零向量")
    axis = axis / norm

    # 3. 将角度转换为弧度
    theta = np.radians(angle_degrees)

    # 4. 提取分量
    ux, uy, uz = axis
    c = np.cos(theta)
    s = np.sin(theta)
    t = 1 - c

    # 5. 构建旋转矩阵
    # 罗德里格公式展开形式
    R = np.array([
        [t * ux ** 2 + c, t * ux * uy - s * uz, t * ux * uz + s * uy],
        [t * ux * uy + s * uz, t * uy ** 2 + c, t * uy * uz - s * ux],
        [t * ux * uz - s * uy, t * uy * uz + s * ux, t * uz ** 2 + c]
    ])

    return R


def main():
    print("--- 3D 旋转矩阵计算器 ---")

    # 获取用户输入
    try:
        axis_input = input("请输入旋转轴向量 (以空格分隔，例如 '0 1 0'): ")
        angle_input = input("请输入旋转角度 (度): ")

        # 解析输入
        axis = [float(x) for x in axis_input.split()]
        angle = float(angle_input)

        if len(axis) != 3:
            print("错误: 旋转轴必须包含 3 个数值 (x y z)")
            return

        # 计算矩阵
        R = get_rotation_matrix(axis, angle)

        # 输出结果
        print("\n计算结果 (保留4位小数):")
        with np.printoptions(precision=4, suppress=True):
            print(R)

        # 验证: 如果输入的是 Y轴 (0 1 0)，可以对比之前的公式
        if axis == [0.0, 1.0, 0.0] or axis == [0, 1, 0]:
            print("\n(注: 这与之前提到的绕Y轴旋转矩阵公式一致)")

    except ValueError:
        print("输入格式错误，请输入数字。")


if __name__ == "__main__":
    main()
