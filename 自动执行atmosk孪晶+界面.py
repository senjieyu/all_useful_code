import subprocess
import os
import sys

# ================= 配置区域 =================
# 如果你的 atomsk 命令不是 'atomsk' (例如是完整路径)，请在这里修改
ATOMSK_CMD = "atomsk"

# 定义文件名，方便统一管理
FILE_CELL = "Cu_cell.xsf"
FILE_MIRROR = "Cu_mirror.xsf"
FILE_FINAL_Y = "Cu_final.cfg"
FILE_ROTATED = "output.cfg"
FILE_MERGED_TEMP = "merged_temp.cfg"
FILE_RESULT = "cu.cfg"


# ================= 工具函数 =================
def run_command(args):
    """执行 Atomsk 命令并处理输出"""
    cmd_str = " ".join(args)
    print(f"\n[执行命令] {cmd_str}")

    try:
        # check=True 会在命令执行出错时抛出异常
        # capture_output=False 让 atomsk 的输出直接显示在终端里
        subprocess.run(args, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"!!! 命令执行出错: {cmd_str}")
        sys.exit(1)


def remove_if_exists(filename):
    """如果文件存在则删除，防止 Atomsk 询问是否覆盖导致脚本卡住"""
    if os.path.exists(filename):
        print(f"[清理旧文件] 删除 {filename}")
        os.remove(filename)


# ================= 主流程 =================
def main():
    print(">>> 开始自动化构建流程...")

    # ---------------------------------------------------------
    # 步骤 1: 创建基础晶胞并扩胞 (3x1x5)
    # ---------------------------------------------------------
    # 命令: atomsk --create fcc 3.615 Cu orient [11-2] [111] [-110] -duplicate 3 1 5 Cu_cell.xsf
    remove_if_exists(FILE_CELL)
    run_command([
        ATOMSK_CMD,
        "--create", "fcc", "3.615", "Cu",
        "orient", "[11-2]", "[111]", "[-110]",
        "-duplicate", "3", "1", "5",
        FILE_CELL
    ])

    # ---------------------------------------------------------
    # 步骤 2: 沿 Y 轴镜像并 Wrap
    # ---------------------------------------------------------
    # 命令: atomsk Cu_cell.xsf -mirror 0 Y -wrap Cu_mirror.xsf
    remove_if_exists(FILE_MIRROR)
    run_command([
        ATOMSK_CMD,
        FILE_CELL,
        "-mirror", "0", "Y",
        "-wrap",
        FILE_MIRROR
    ])

    # ---------------------------------------------------------
    # 步骤 3: 第一次拼接 (沿 Y 轴堆叠原始和镜像)
    # ---------------------------------------------------------
    # 命令: atomsk --merge Y 2 Cu_cell.xsf Cu_mirror.xsf Cu_final.cfg
    remove_if_exists(FILE_FINAL_Y)
    run_command([
        ATOMSK_CMD,
        "--merge", "Y", "2",
        FILE_CELL, FILE_MIRROR,
        FILE_FINAL_Y
    ])

    # ---------------------------------------------------------
    # 步骤 4: 旋转操作 (沿 Z 轴旋转 180 度)
    # ---------------------------------------------------------
    # 命令: atomsk Cu_final.cfg -rotate Z 180 output.cfg
    remove_if_exists(FILE_ROTATED)
    run_command([
        ATOMSK_CMD,
        FILE_FINAL_Y,
        "-rotate", "Z", "180",
        FILE_ROTATED
    ])

    # ---------------------------------------------------------
    # 步骤 5: 第二次拼接 (沿 Z 轴堆叠未旋转和旋转后的部分)
    # ---------------------------------------------------------
    # 命令: atomsk --merge Z 2 Cu_final.cfg output.cfg merged_temp.cfg
    remove_if_exists(FILE_MERGED_TEMP)
    run_command([
        ATOMSK_CMD,
        "--merge", "Z", "2",
        FILE_FINAL_Y, FILE_ROTATED,
        FILE_MERGED_TEMP
    ])

    # ---------------------------------------------------------
    # 步骤 6: 最终 Wrap 处理 (标准化坐标)
    # ---------------------------------------------------------
    # 命令: atomsk merged_temp.cfg -wrap cu.cfg
    remove_if_exists(FILE_RESULT)
    run_command([
        ATOMSK_CMD,
        FILE_MERGED_TEMP,
        "-wrap",
        FILE_RESULT
    ])

    print(f"\n>>> 流程结束！最终文件已生成: {FILE_RESULT}")

    # ---------------------------------------------------------
    # 可选: 清理中间文件
    # ---------------------------------------------------------
    # 如果你想保留中间文件，请注释掉下面这几行
    clean = input("\n是否删除中间过程文件 (y/n)? ")
    if clean.lower() == 'y':
        intermediate_files = [FILE_CELL, FILE_MIRROR, FILE_FINAL_Y, FILE_ROTATED, FILE_MERGED_TEMP]
        for f in intermediate_files:
            remove_if_exists(f)
        print(">>> 中间文件已清理。")


if __name__ == "__main__":
    main()
