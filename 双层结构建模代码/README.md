# 双层/梯度纳米孪晶铜建模代码

本仓库用于生成单层晶粒、纳米孪晶块、自定义多层叠层，以及 AB、ABAB、TTG、GSG、DG 和纯相参照模型。输出包含 Atomsk 可读的 CFG/XSG 结构、LAMMPS `wrapped.lmp` 数据文件及完整元数据。

代码采用**一个配置文件驱动**的方式。日后需要改变晶格常数、盒子大小、层厚、孪晶间距、晶粒模板、模型清单、输出目录或默认任务，只改 `config.toml`；不必修改 Python 主程序。

## 文件说明

| 文件 | 用途 |
|---|---|
| `dual_gradient_cu_builder.py` | 主程序：生成块、合并层并批量构建模型。 |
| `config.toml` | 唯一参数入口，所有参数均有中文注释。 |
| `run_builder.bat` | Windows 双击运行入口，按 `config.toml` 的 `[run]` 执行。 |
| `test_dual_gradient_cu_builder.py` | 算法单元测试。 |
| `output/` | 默认输出目录；首次建模时自动生成，不随代码发布。 |

## 环境准备

1. 安装 Python 3.8 或更高版本。
2. 安装 [Atomsk](https://atomsk.univ-lille.fr/)。
3. 在 `config.toml` 的 `[paths]` 中设置 Atomsk：
   - 已将 Atomsk 加入系统 `PATH`：`atomsk = "auto"`；
   - 未加入 `PATH`：填写完整路径，例如 `atomsk = "D:/software/Atomsk/atomsk.exe"`。
4. 在本目录打开 PowerShell 或 CMD。

程序不依赖第三方 Python 包。Python 3.11 及以上直接使用内置 TOML 读取器；Python 3.8–3.10 自动使用脚本内的兼容读取器。

## 最快使用方式

1. 打开 `config.toml`。
2. 在 `[paths]` 中确认 `atomsk` 和 `output_root`。
3. 在 `[build]` 中确认晶格常数 `a0_angstrom = 3.633`、面内尺寸和种子参数。
4. 在 `[run]` 中设置本次任务。
5. 双击 `run_builder.bat`，或在终端运行：

```powershell
python dual_gradient_cu_builder.py
```

先只检查配置是否被正确读取：

```powershell
python dual_gradient_cu_builder.py --show-config
```

## 配置文件说明

### 1. 全局几何与输出

- `[paths].output_root`：模型输出根目录。可使用相对路径，例如 `output`。
- `[build].a0_angstrom`：FCC Cu 晶格常数，当前固定为 `3.633` Å。
- `[build].box_xy_nm`：面内盒子尺寸。增大该值会增加原子数、内存和建模时间。
- `[build].twin_nx`、`twin_ny`：孪晶种子面内重复数。

### 2. 孪晶、层厚与晶粒模板

- `[twin_n111]`：`T4`、`T6` 等代码对应一个孪晶周期的 `{111}` 原子层数；脚本根据 `a0` 自动计算实际间距。
- `[layer_thickness_nm]`：AB、ABAB 和纯相基线的层厚度。
- `[grain_templates.*]`：`nx × ny` 是面内 Voronoi 晶粒排布。例如 `G3x2` 代表 6 个晶粒。

修改或新增晶粒模板后，`plan_models` 和 `[run]` 中引用的名称必须一致。修改或新增 T 编号后，同样应使用一致的孪晶代码。

### 3. 批量模型清单

每个 `[[plan_models]]` 段定义一个双层模型：

```toml
[[plan_models]]
name = "TTG-2"
family = "TTG"
phase = "P0"
grain_a = "G3x2"
grain_b = "G3x2"
twin_a = "T10"
twin_b = "T20"
```

可直接增加、删除或修改这些段落。`family` 用于输出分类和筛选，`phase` 用于分批建模，`grain_a/B` 与 `twin_a/B` 分别定义两种组成层。

## 三种运行模式

在 `config.toml` 的 `[run]` 段设置 `command`。

### A. 批量建模：`command = "plan"`

默认模式。用于按 `[[plan_models]]` 批量生成 AB、ABAB 与纯相参照结构。

```toml
[run]
command = "plan"
plan_phase = "P0"       # all、P0、P1 或 P2
include_pure = true
families = ["TTG", "DG"]
models = []
variants = ["AB"]       # [] 为同时构建 AB 与 ABAB
force = false
```

- `families`、`models`、`variants` 为 `[]` 时不筛选。
- 将 `force = true` 用于参数改动后强制重建；否则已有结构会被保留。

### B. 单层块：`command = "block"`

用于先检查某个晶粒/孪晶组合：

```toml
[run]
command = "block"
block_grain = "G3x2"
block_twin = "T15"       # 无孪晶时写 "none"
block_thickness_nm = 24.0
```

### C. 自定义叠层：`command = "merge"`

用于任意层数和任意组合：

```toml
[run]
command = "merge"
merge_name = "custom_AB"
merge_family = "CUSTOM"
merge_phase = "custom"
merge_layers = ["G3x2:T10:24", "G3x2:T20:20"]
```

每个 `merge_layers` 条目的格式为 `晶粒模板:孪晶代码:厚度_nm`。列表顺序即 Atomsk 的 Z 向合并顺序。

## 命令行临时覆盖

日常参数应在 `config.toml` 修改。临时试算时，也可在不改配置文件的情况下使用命令行子命令：

```powershell
# 使用当前配置中的几何参数，临时构建一个块
python dual_gradient_cu_builder.py block --grain G3x2 --twin T15 --thickness 24

# 临时构建一组 P0 的 TTG-2 的 AB 模型
python dual_gradient_cu_builder.py plan --phase P0 --models TTG-2 --variants AB

# 使用另一份配置
python dual_gradient_cu_builder.py --config my_config.toml
```

命令行中给出的参数只影响这一次运行；配置文件不会被改写。

## 输出结构

以 `output_root = "output"` 为例：

```text
output/
├─ cache/blocks/             # 单层块缓存，避免重复构建
├─ models/<family>/<name>/<variant>/
│  ├─ wrapped.lmp            # LAMMPS 数据文件
│  ├─ wrapped.cfg            # Atomsk CFG 文件
│  ├─ merged.cfg             # 合并后的中间结构
│  └─ metadata.json          # 结构参数与来源记录
└─ manifests/
   ├─ plan_reference.json    # 本次参数、孪晶间距和模型清单
   └─ build_manifest_*.json  # 批量构建结果
```

推荐将 `wrapped.lmp` 用作 LAMMPS 后续弛豫和拉伸的初始构型，并保留同目录的 `metadata.json` 与 `plan_reference.json` 以追溯参数。

## 验证与常见问题

运行算法测试：

```powershell
python test_dual_gradient_cu_builder.py
```

- **找不到 Atomsk**：检查 `config.toml` 的 `atomsk` 路径，或将 Atomsk 加入 `PATH` 后写 `auto`。
- **改完参数后结果未变化**：将 `[run].force` 改为 `true`，或修改 `output_root` 使用新目录。
- **模型清单报模板/孪晶代码错误**：确认 `[[plan_models]]` 的 `grain_a/B` 和 `twin_a/B` 已在 `[grain_templates]` 与 `[twin_n111]` 中定义。
- **内存不足**：先减小 `box_xy_nm`，或只构建 `variants = ["AB"]`、单个 `models` 进行检查。

## 发布建议

发布时保留本目录下的 5 个文件：`dual_gradient_cu_builder.py`、`config.toml`、`run_builder.bat`、`test_dual_gradient_cu_builder.py` 和本 README。不要发布 `output/` 中的大型模型文件；需要共享某个已生成结构时，单独打包对应模型目录及其 `metadata.json`。
