# 推文代码：Atomsk + LAMMPS 晶界工作流

本目录整理自两篇微信推文正文中可见的 Atomsk 与 LAMMPS 命令，用于构建
FCC Cu 的 Σ17{223}[110] 对称倾侧晶界，并以刚体平移筛选低晶界能构型。

## 入口

- Windows：`scripts/generate_sigma17_stgb.bat`
- Linux/macOS：`scripts/generate_sigma17_stgb.sh`
- LAMMPS 弛豫：`lammps/in.optimize_stgb.lmp`

运行前需自行安装 Atomsk、LAMMPS，并向 `lammps/` 放入与 `pair_style eam`
兼容的铜势文件（示例名：`Cu_u3.eam`）。

## 来源

- https://mp.weixin.qq.com/s/u6qvvsnszPU6pr8u4O0Tgw
- https://mp.weixin.qq.com/s/ZXrnbjKbsKerm_ZVDSoDFQ

推文网盘附件未被复制；此处仅包含从正文可见命令重新排版和参数化的代码。
