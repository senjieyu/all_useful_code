# -*- coding: utf-8 -*-
"""
Atomsk + Python 生成 BCC Fe Σ3 {112} 孪晶结构（输出 POSCAR/.vasp）
兼容 Atomsk 0.13.1

输出：
- Fe_BCC_S3_112_twin.vasp
"""

import os
import shutil
import subprocess


def run(cmd):
    print("  执行:", " ".join(cmd))
    p = subprocess.run(cmd, text=True, capture_output=True)
    if p.stdout:
        print(p.stdout.rstrip())
    if p.stderr:
        print(p.stderr.rstrip())
    if p.returncode != 0:
        raise RuntimeError(f"命令失败: {' '.join(cmd)}")


def main():
    atomsk = shutil.which("atomsk") or shutil.which("atomsk.exe")
    if atomsk is None:
        raise RuntimeError("找不到 atomsk/atomsk.exe，请把 Atomsk 加入 PATH 或放到脚本目录。")

    print("=" * 60)
    print("BCC Fe Σ3 {112} twin（mirror + merge） -> POSCAR/.vasp")
    print("=" * 60)

    # 检查 atomsk
    run([atomsk, "--version"])

    # ---------------------------------------------------------
    # 1) 创建 BCC Fe
    #    关键：Y = [112] 作为(112)法向
    #    Atomsk 0.13.1 对这组基矢要求较严：Z 用 [22-2]（你实际跑过可行）
    # ---------------------------------------------------------
    run([
        atomsk,
        "--create", "bcc", "2.87", "Fe",
        "orient", "[-110]", "[112]", "[22-2]",
        "-duplicate", "1", "14", "1",
        "Fe_A.xsf"
    ])

    # 2) 镜像生成孪晶晶粒
    run([atomsk, "Fe_A.xsf", "-mirror", "0", "Y", "-wrap", "Fe_B.xsf"])

    # 3) 沿 Y 合并形成孪晶界面
    run([atomsk, "--merge", "Y", "2", "Fe_A.xsf", "Fe_B.xsf", "Fe_twin.xsf"])

    # 4) 输出 POSCAR
    run([atomsk, "Fe_twin.xsf", "pos"])

    # 5) 重命名
    out_name = "Fe_BCC_S3_112_twin.vasp"
    if os.path.exists(out_name):
        os.remove(out_name)
    os.rename("POSCAR", out_name)

    print("=" * 60)
    print("完成输出：", out_name)
    print("=" * 60)


if __name__ == "__main__":
    main()
