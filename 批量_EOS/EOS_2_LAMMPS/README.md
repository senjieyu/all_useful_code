# EOS 批量计算作业流程使用说明

本说明文档介绍了如何使用批量处理脚本将 `.xyz` 文件中的多个结构自动生成独立的 EOS 计算任务，并批量提交到计算集群。

## 1. 准备工作

请确保以下文件都在同一个目录下：

*   **输入结构文件**: `CU331.xyz` (或您的其他包含多帧结构的 .xyz 文件)
*   **必要配置文件**:
    *   `mlip.ini`
    *   `relax_isif2.in`
    *   `current.mtp`
*   **作业脚本模板**: `submit_eos_isif2.sh`
*   **Python 脚本**:
    *   `rand_dis_battery_data.py` (核心数据生成脚本)
    *   `setup_batch.py` (用于初始化文件夹和配置)
    *   `submit_batch.py` (用于批量提交作业)

## 2. 运行步骤

### 第一步：生成作业环境

使用 `setup_batch.py` 脚本来处理输入文件，它会自动为每一帧结构创建独立的文件夹，并配置好所需的所有文件。

**命令格式：**
```bash
python setup_batch.py <输入xyz文件> <输出总文件夹名>
```

**示例：**
```bash
python setup_batch.py CU331.xyz batch_work
```

**执行结果：**
*   在当前目录下创建一个名为 `batch_work` 的总文件夹。
*   在 `batch_work` 中，为每个结构生成子文件夹（如 `struct_0`, `struct_1`...）。
*   每个子文件夹内包含：
    *   `eos_data/`: 包含不同缩放比例的 LAMMPS data 文件。
    *   `submit.sh`: 针对该任务定制的提交脚本。
    *   必要的配置文件 (`mlip.ini`, `relax_isif2.in`, `current.mtp`)。
*   在 `batch_work/submission_scripts/` 下生成所有任务的启动脚本。

---

### 第二步：批量提交作业

使用 `submit_batch.py` 脚本一次性提交所有生成的作业。

**命令格式：**
```bash
python submit_batch.py <提交脚本文件夹路径>
```

**示例：**
```bash
python submit_batch.py batch_work/submission_scripts
```

**执行结果：**
*   脚本会遍历指定文件夹中的所有 `.sh` 文件。
*   自动调用 `bsub` 命令将任务提交到队列。
*   终端会显示提交成功或失败的统计信息。

## 3. 结果查看

作业提交后，计算结果将保存在各个结构的独立文件夹中。

例如，对于第 0 帧结构：
*   路径：`batch_work/struct_0/`
*   日志文件：`batch_work/struct_0/eos_logs/`
*   结果数据：`batch_work/struct_0/eos_results/`

## 4. 常见问题

*   **找不到文件错误**：请确保运行脚本时，所有“准备工作”中列出的文件都在当前目录下。
*   **权限错误**：如果提示 `bash` 命令未找到或权限不足，请确保在 Linux 环境或支持 bash 的环境（如 WSL）中运行提交脚本。
