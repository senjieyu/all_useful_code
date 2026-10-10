# Random Forest Diagnostic Suite

将微信推文中可见的随机森林二分类分析思路整理为可运行脚本：从训练、调参、解释性到业务阈值决策，一次输出 18 类图表。

## 原文图文

- 原文：[机器学习——随机森林+PDP + ICE+ROC +Shap+代价曲线+决策曲线+KS曲线+Lift+六视角横评模型](https://mp.weixin.qq.com/s/Lg4tk91s4eVcvmyw8dW6fg)
- 本目录仅保存重新组织的可运行实现、使用说明与原文链接；不复制推文完整图文、图片或需要私信获取的数据。

原文涵盖混淆矩阵、OOB/CV 调参、PDP/ICE、学习曲线、ROC/PR、校准、SHAP、成本曲线、决策曲线、KS、Lift 和误分类诊断。本实现保留这些分析目标，并修正阈值成本计算为：漏判正类（FN）与误判负类（FP）的显式成本。

## 安装

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
python -m pip install -r requirements.txt
```

## 快速演示

不需要外部数据，先生成一个合成的、信贷风险风格二分类数据集：

```bash
python random_forest_diagnostic_suite.py --demo --output demo_results
```

若没有安装 `shap`，可先跳过 SHAP 图：

```bash
python random_forest_diagnostic_suite.py --demo --skip-shap
```

## 使用自己的 CSV

要求目标列是二分类；非 `0/1` 的两个类别会自动编码。除目标列外仅使用数值特征。

```bash
python random_forest_diagnostic_suite.py \
  --data your_data.csv \
  --target credit \
  --output rf_results
```

输出目录包含 `01` 至 `18` 编号的 PNG。阈值的选择严格使用训练集 OOF 概率，测试集只用于最终报告，避免测试集泄漏。业务上线前仍须根据实际 FN/FP 成本、类别定义、时间切分和校准结果重新确认参数。
