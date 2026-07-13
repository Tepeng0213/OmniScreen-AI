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

## 蛋白/多肽线典型产出

| 文件 | 阶段 | 说明 |
|------|------|------|
| `mutation_scores.csv` | Module 2 | ESM-2 CDR 突变 ΔLL 打分 |
| `ppi_docking_scores.csv` | Module 3 | 蛋白-蛋白界面评分 |
| `pe_docking/*.pdb` | Module 3 | 纳米抗体结构（extended-CA 或 ESMFold） |
| `af3_pe_metrics.csv` | Module 4 | AF3 纳米抗体–PD-L1 ipTM / pTM |
| `af3_pe_complexes/*.cif` | Module 4 | AF3 最佳复合物结构 |
| `figures/fig5*.png`, `fig6*.png`, `fig_pe_af3_*` | Module 4/6 | CDR 热图、PPI、AF3 排名/3D |
| `af3_server/` | Module 4 准备 | AlphaFold Server 上传 JSON（PE/NA/SM 材料） |

详见 [PE 模块说明](../../docs/workflows/PE_MODULES.md)。

## 核酸线典型产出

| 文件 | 阶段 | 说明 |
|------|------|------|
| `CD274_mRNA.fasta` | Module 1 | PD-L1 mRNA 序列 |
| `sirna_candidates.csv` | Module 2 | siRNA 候选 + ViennaRNA 打分 |
| `offtarget_filtered.csv` | Module 3 | 脱靶过滤后候选 |
| `sirna_offtarget.sam` | Module 3 | Bowtie2 比对结果（chr22 演示） |
| `af3_na_metrics.csv` | Module 4 | AF3 蛋白–siRNA ipTM / pTM |
| `af3_na_complexes/*.cif` | Module 4 | AF3 最佳复合物结构 |
| `thermodynamics.csv` | Module 5 | 杂交 ΔG / 发夹 ΔG |
| `figures/fig7*.png`, `fig8*.png`, `fig_na_*.png` | Module 4–6 | 效能、脱靶、AF3、热力学图 |

详见 [NA 模块说明](../../docs/workflows/NA_MODULES.md)。

## 跨平台数据交接

Colab Module 1–3 完成后，结果已 push 到 GitHub。RunPod 上 `git clone` 同一仓库即可继续 Module 4–5。

## 蛋白质/多肽线典型产出

| 文件 | 阶段 | 说明 |
|------|------|------|
| `af3_pe_metrics.csv` / `af3_pe_complexes/` | Module 4 | AF3 复合物与 ipTM |
| `ppi_mmgbsa_summary.csv` | Module 5 | MM-GBSA 结合自由能排序 |
| `ppi_energy_decomposition.csv` | Module 5 | 残基级 VDW/ELE 拆解 |
| `figures/fig_pe_mmgbsa_*.png` | Module 5 | ΔG 排序 / 残基贡献 / vs ipTM |
