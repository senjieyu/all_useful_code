# VASP Elastic Constants Batch Calculation

这是一个用于批量使用 VASP 计算 `.xyz` 结构文件弹性常数的工具包。
流程包括：将 XYZ 转为 POSCAR -> 生成独立作业文件夹 -> 批量提交 -> 收集 OUTCAR 中的弹性模量信息。

## 1. 文件夹内容

*   `setup_vasp_elastic.py`: 初始化脚本。
*   `submit_batch.py`: 批量提交脚本（通用）。
*   `collect_vasp_elastic.py`: 结果收集脚本。
*   `INCAR_relax`: 结构弛豫参数。
*   `INCAR_elastic`: 弹性计算参数 (IBRION=6)。

## 2. 使用步骤

### 第一步：生成作业

```bash
# 用法: python setup_vasp_elastic.py <输入xyz> <输出目录>
python setup_vasp_elastic.py ../CU331.xyz Vasp_Work
```

**注意**：脚本会在每个子文件夹中生成 `submit.sh`，但**不会**复制 `POTCAR`。
**必须**在提交前将 `POTCAR` 放入每个文件夹，或者修改 `setup_vasp_elastic.py` 让其自动复制。

批量复制 POTCAR 示例：
```bash
for d in Vasp_Work/struct_*; do cp /path/to/your/POTCAR $d/; done
```

### 第二步：提交作业

```bash
python submit_batch.py Vasp_Work/submission_scripts
```
作业会自动分两步运行：
1.  结构弛豫 (`INCAR_relax`) -> 更新 `POSCAR`
2.  弹性计算 (`INCAR_elastic`) -> 输出 `OUTCAR`

### 第三步：收集结果

```bash
python collect_vasp_elastic.py Vasp_Work vasp_results
```
结果将从 `OUTCAR` 中提取 "TOTAL ELASTIC MODULI" 部分并保存在 `vasp_results` 文件夹中。

## 3. INCAR 配置说明

您可以直接编辑 `INCAR_relax` 和 `INCAR_elastic` 来修改计算参数（如 ENCUT, KSPACING 等）。
默认配置：
*   Relax: IBRION=2, ISIF=3
*   Elastic: IBRION=6, ISIF=3, NFREE=2
