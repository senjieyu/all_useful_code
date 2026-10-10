# Σ17{223}[110] 铜晶界：Atomsk、LAMMPS 与 aimsgb

本目录包含两套可配合使用的晶界建模代码：

- `wechat_article_code/`：根据微信推文正文整理的 Atomsk + LAMMPS 工作流。它建立 FCC Cu 的 Σ17{223}[110] 对称倾侧晶界，通过刚体平移和能量最小化筛选低能构型。
- `github_code/`：aimsgb 上游源码快照与基于 JSON 配置的一键晶界生成器。

## 已验证的环境与范围

以下流程已在 Windows 上使用 Atomsk 0.13.1 和 LAMMPS 7 Aug 2019 完成实际冒烟测试：Atomsk 生成 27,200 原子双晶；LAMMPS 读取两种 atom type、将它们共同映射至 `Cu_u3.eam`，并完成缩小模型的 2×2 刚体平移、重叠删除和能量最小化。

几何关系为：旋转轴 `[110]`、`cos(theta)=1/17`（`theta≈86.63°`），合并方向对应 `{223}` 晶界面。Atomsk 将第二晶粒临时标记为 Al 以保留 atom type 2；LAMMPS 输入会将两个 type 的质量都重设为 Cu。

## 方式一：Atomsk + LAMMPS（Windows）

1. 安装 Atomsk 和 LAMMPS，并编辑 `wechat_article_code/scripts/generate_sigma17_stgb.bat` 中的 `ATOMSK_BIN`，使其指向本机 `atomsk.exe`。
2. 在 `wechat_article_code/scripts/` 目录运行：

   ```bat
   generate_sigma17_stgb.bat
   ```

   它会生成 `sigma17_stgb_initial.lmp`。
3. 将该数据文件复制到 `wechat_article_code/lammps/`，并将与 LAMMPS `pair_style eam` 兼容的 Cu 势文件复制到同一目录、命名为 `Cu_u3.eam`。
4. 在 `wechat_article_code/lammps/` 运行：

   ```bat
   lmp_serial.exe -in in.cohesive_energy.lmp
   lmp_serial.exe -in in.optimize_stgb.lmp
   ```

   `in.cohesive_energy.lmp` 还需要一个 FCC Cu 单晶数据文件 `bulk_cu.lmp`。将该计算得到的内聚能回填到 `in.optimize_stgb.lmp` 的 `cohesive_energy` 后，再运行晶界扫描。
5. 在 Linux、WSL 或 Git Bash 中按能量排序：

   ```bash
   ../scripts/sort_gb_energies.sh log.lammps > gb_energy_sorted.txt
   ```

输出的 `result/after_minimize_<i>_<j>.lmp` 是各刚体平移构型弛豫后的候选结构；选择最低晶界能对应的文件用于后续模拟。

## 方式二：aimsgb 配置驱动生成器

进入 `github_code/csl_grain_boundary_builder/`，使用 Python 3.9+：

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install -e .
python -m csl_gb_builder.cli build configs/sigma17_110_223_cu.json --dry-run
python -m csl_gb_builder.cli build configs/sigma17_110_223_cu.json
```

可修改 `configs/sigma17_110_223_cu.json` 中的元素、晶格常数、Σ 值、旋转轴、晶界面、超胞大小、真空层和输出格式。图形入口为 `run_builder.py`；安装 PyInstaller 后，在 Windows 双击 `build_exe.bat` 可打包 EXE，配置文件保持在 EXE 旁边，后续修改参数不必重新打包。

## 来源与许可证

- 微信推文：[任意CSL值晶界建模(一)](https://mp.weixin.qq.com/s/u6qvvsnszPU6pr8u4O0Tgw)、[任意CSL值晶界建模(二)](https://mp.weixin.qq.com/s/ZXrnbjKbsKerm_ZVDSoDFQ)
- aimsgb 项目：[GitHub](https://github.com/ksyang2013/aimsgb)，[网页生成器](https://aimsgb.org/)

`github_code/aimsgb_source/` 保留 aimsgb 1.1.1 上游快照和版权声明。仅可按其教育、科研、非营利许可使用、修改和分发；商业使用须另行取得授权。详见 `github_code/THIRD_PARTY_NOTICES.md`。
