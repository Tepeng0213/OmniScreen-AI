# 初始候选库

本目录存放海选阶段的原始数据（不纳入 Git 追踪）。

## 小分子线

- `initial_compounds.smi` — 空格分隔：`SMILES  MOL_ID`
- 来源：ChEMBL / ZINC20 子集，或 notebook 内置测试集

## 蛋白/多肽线

- `pd_l1_nanobody_seed.fasta` — 抗体种子序列（如 KN035）
- 突变库由 notebook 内生成

## 核酸线

- `CD274_mRNA.fasta` — 人类 CD274 (PD-L1) mRNA（NCBI NM_014143）
- siRNA 候选由 notebook 滑动窗口自动生成
