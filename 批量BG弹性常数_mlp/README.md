# Elastic Constants Batch Calculation

这是一个独立的工具包，用于批量计算 `.xyz` 结构文件中的弹性常数。
它基于 LAMMPS 和 MLIP 势函数，自动进行结构弛豫、施加应变并计算 Stiffness Tensor (Cij) 及相关模量（Bulk, Shear, Poisson's Ratio 等）。

## 1. 文件夹内容

*   `setup_elastic.py`: 初始化脚本，将 xyz 文件拆分为独立作业。
*   `submit_batch.py`: 批量提交作业脚本。
*   `collect_elastic.py`: 结果收集脚本。
*   `calc_elastic_single.py`: 核心计算逻辑（在每个作业中运行）。
*   `mlip.ini` & `current.mtp`: 势函数文件。

## 2. 环境要求

*   Python 3 (需要安装 `ase`, `numpy`)
*   LAMMPS (需要编译 MLIP 接口)
*   运行环境支持 `bsub` (LSF 作业调度系统)

## 3. 使用步骤

假设你有一个结构文件 `CU331.xyz`（可以位于上一级目录或任意位置）。

### 第一步：生成作业

运行 `setup_elastic.py`。这会创建一个工作目录，为每个结构生成独立的计算文件夹。

```bash
# 用法: python setup_elastic.py <输入xyz文件> <输出工作目录>
python setup_elastic.py ../CU331.xyz Elastic_Work
```

*   `Elastic_Work`: 你自定义的文件夹名，所有的计算任务都会放在这里。

### 第二步：提交作业

运行 `submit_batch.py`。它会遍历工作目录下的 `submission_scripts` 并提交所有任务。

```bash
# 用法: python submit_batch.py <脚本目录>
python submit_batch.py Elastic_Work/submission_scripts
```

### 第三步：收集结果

等待作业运行完成后，运行 `collect_elastic.py`。它会将每个任务的计算结果汇总到一个文件夹中。

```bash
# 用法: python collect_elastic.py <工作目录> <结果输出目录>
python collect_elastic.py Elastic_Work elastic_results
```

## 4. 结果说明

结果文件夹 (`elastic_results`) 中将包含多个 `.txt` 文件（如 `struct_0.txt`, `struct_1.txt`...）。

每个文件内容示例：
```text
Stiffness Tensor Cij (GPa)
    180.00     127.00     127.00       0.00       0.00       0.00
    127.00     180.00     127.00       0.00       0.00       0.00
    ...

Bulk Modulus Kv:   145.00 GPa
Shear Modulus Gv:   57.00 GPa
Elastic Anisotropy:  1.51
Poisson's Ratio:     0.34
```
