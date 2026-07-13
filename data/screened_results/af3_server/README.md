# AlphaFold 3 Server 提交材料

本目录由 OmniScreen CPU 阶段（Module 0–3）结果自动整理，供 [AlphaFold Server](https://alphafoldserver.com) 上传使用。

> **重要**：Server **没有**公开 API，不能从 Cursor/Colab 直接遥控。流程是：notebook/本目录生成 JSON → 浏览器 Upload JSON → 下载 zip → 再回 notebook 做可视化。

## 目录结构

```
af3_server/
├── README.md                 # 本说明
├── pe/                       # 蛋白/多肽：纳米抗体–PD-L1
│   ├── batch_top5.json       # ★ 推荐：一次上传 5 个 job
│   ├── PE_*.json             # 单 job 备份
│   ├── manifest.csv          # 候选元数据
│   └── pdl1_4zqk_chainB.fasta
├── na/                       # 核酸：PD-L1 + siRNA
│   ├── batch_top5.json       # ★ 推荐：蛋白 + antisense + sense（双链）
│   ├── batch_top5_antisense_only.json  # 更轻：仅 antisense
│   ├── NA_*.json
│   └── manifest.csv
└── sm/                       # 小分子：材料包（Server 能力有限）
    ├── SM_PDL1_5N2F_dimer_baseline.json  # PD-L1 二聚体基线（无自定义配体）
    ├── manifest_top5_ligands.csv
    ├── top5_ligands.smi      # 供 OpenMM MD / 外部对接
    └── pdl1_5n2f_chains.fasta
```

## 上传步骤（PE / NA）

1. 打开 https://alphafoldserver.com 并登录 Google 账号  
2. 点击 **Upload JSON**  
3. 上传：
   - PE：`pe/batch_top5.json`
   - NA：`na/batch_top5.json`（若配额紧，用 `batch_top5_antisense_only.json`）
4. 在 History 中检查草稿 → **Run**（消耗每日免费额度）  
5. 完成后下载 zip，解压到例如：

```text
data/screened_results/af3_results/pe/
data/screened_results/af3_results/na/
```

zip 内通常含：
- `*_model.cif` / PDB
- `*_confidences.json`（含 ipTM / pTM / PAE）
- `*_job_request.json`（复现用输入）

## 各线路说明

### PE（优先，最适合 Server）

| 项目 | 内容 |
|------|------|
| 复合物 | Top 5 纳米抗体突变体 + PD-L1（4ZQK chain B，106 aa） |
| 候选来源 | `ppi_docking_scores.csv` 前 5：`CDR3_P98KV/A/L/I/P` |
| 预期产出 | 抗体–抗原界面结构 + ipTM / PAE |

### NA（适合 Server）

| 项目 | 内容 |
|------|------|
| 复合物 | PD-L1 + siRNA antisense（+ sense 双链） |
| 候选来源 | `offtarget_filtered.csv` 中 `passed_offtarget` Top 5 |
| 注意 | 生物学上 siRNA 主要靶向 **mRNA**；此处蛋白–核酸 job 用于结构/空间验证演示，解读时需谨慎 |

### SM（Server 几乎跑不了自定义配体）

AlphaFold Server **只允许白名单 CCD 配体**（ATP、HEM、NAD 等），**不接受任意 SMILES**。  
你们的 Top 配体（`CHEMBL19019_div` 等）无法直接提交。

已准备：

| 文件 | 用途 |
|------|------|
| `SM_PDL1_5N2F_dimer_baseline.json` | 可选：跑 PD-L1 二聚体基线结构 |
| `top5_ligands.smi` + `manifest_top5_ligands.csv` | 交给 **OpenMM MD / MM-PBSA**（SM Module 4–5 正路） |
| `5N2F.pdb`（仓库内） | MD 对接受体 |

**SM 建议路径**：AF3 Server 可跳过自定义配体 job → 用 RunPod 做 OpenMM。

## 配额建议

| 线路 | job 数 | 建议 |
|------|--------|------|
| PE | 5 | 先跑 |
| NA | 5（或 antisense-only 5） | 次日或配额够时再跑 |
| SM | 0–1（仅二聚体基线） | 可选 |

## 跑完后回 Cursor 做什么

1. 把 zip 解压到 `data/screened_results/af3_results/{pe,na}/`  
2. 在 PE/NA notebook Module 4/6 中解析 CIF + confidence → 出 3D / 雷达图  
3. `export_for_local_sync()` 写回本地（若在 Colab 跑可视化）

---

*生成自 OmniScreen CPU 阶段结果 · 2026-07*
