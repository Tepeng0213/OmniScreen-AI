# 筛选结果输出

各 workflow notebook 各模块的输出文件存放于此。

**持久化规则**：CSV / PNG 等小型结果会通过 `persist_to_github()` 自动提交到 [GitHub](https://github.com/Tepeng0213/OmniScreen-AI)。MD 轨迹（`.dcd`）过大，不纳入 Git。

## 小分子线典型产出

| 文件 | 阶段 | 说明 |
|------|------|------|
| `chemical_space_props.csv` | Module 2 | 理化性质 + 过滤标签 |
| `docking_scores.csv` | Module 3 | Vina 对接打分排序 |
| `top10_ligands.sdf` | Module 3 | Top 10 配体（交接 RunPod） |
| `md_rmsd.csv` | Module 4 | MD 轨迹分析 |
| `mmpbsa_results.csv` | Module 5 | 结合自由能 |

## 跨平台数据交接

Colab Module 1–3 完成后，结果已 push 到 GitHub。RunPod 上 `git clone` 同一仓库即可继续 Module 4–5。
