# OmniScreen-AI

**OmniScreen-AI: A Multi-Modal High-Throughput Drug Screening & Dynamic Verification Platform for PD-L1**

面向肿瘤免疫靶点 **PD-L1** 的全模态（小分子 / 大环肽 / 核酸）高通量筛选与动态动力学验证平台。

## 项目亮点

- **三模态并行漏斗**：化学空间（小分子）、序列空间（抗体/多肽）、RNA 空间（siRNA/Aptamer）
- **四段闭环验证**：AI 初筛 → 结构对接 → MD 动力学 → 自由能计算
- **Notebook 驱动**：全部逻辑集中在 3 个 workflow notebook，模块化分块，Colab 一键复现
- **算力分层**：Colab CPU（初筛 + 对接）→ RunPod GPU（MD + 自由能）

## 仓库结构

```
OmniScreen-AI/
├── notebooks/
│   ├── OmniScreen_SM_Workflow.ipynb          # 小分子全流程（ChEMBL 真实化合物库）
│   ├── OmniScreen_PE_Workflow.ipynb   # 蛋白质/多肽全流程
│   └── OmniScreen_NA_Workflow.ipynb   # 核酸药物全流程
├── data/
│   ├── receptor/          # PD-L1 受体 PDB（5N2F, 4ZQK 等）
│   ├── raw_libraries/     # 初始候选库（SMILES / FASTA / mRNA）
│   └── screened_results/  # 各阶段筛选输出
├── docker/                # 可选：标准化运行环境
├── pyproject.toml
└── requirements.txt
```

## 快速开始（Cursor + Colab）

1. Clone 本仓库，用 Cursor 打开
2. 安装 [Google Colab](https://marketplace.visualstudio.com/items?itemName=Google.colab) 与 [Notebook MCP](https://marketplace.visualstudio.com/items?itemName=olavovieiradecarvalho.notebook-mcp-server) 扩展
3. **配置 GitHub Token**（推送数据必需）：
   - 在 [GitHub Settings → Tokens](https://github.com/settings/tokens) 创建 PAT（勾选 `repo` 权限）
   - **Cursor + Colab 扩展**：运行 notebook 中「配置 GitHub Token」cell，手动粘贴 token（Colab Secrets 仅在网页端有效，Cursor 内核读不到）
   - **Colab 网页端**：可在 🔑 Secrets 添加 `GITHUB_TOKEN`，代码会自动读取
4. 打开 `notebooks/OmniScreen_SM_Workflow.ipynb`，选择 **Colab** 内核
5. 按模块顺序运行 cell（Module 0 → Module 6）

> **算力建议**：Module 1–3 在 Colab CPU/GPU 免费层即可；Module 4–5（OpenMM / 自由能）建议在 RunPod 等 GPU 平台运行，详见各 notebook 内说明。

## 数据持久化策略（Colab + 本地自动同步）

**Colab 在云端运行，无法直接写本地磁盘。** 采用「云端计算 → 打包输出 → Agent 写回本地」：

```
Colab 运行模块
  → export_for_local_sync() 将 data/ 文件 base64 输出到 cell
  → Cursor Agent 解析输出标记
  → 自动写入本地 /Users/schmit/Documents/OmniScreen-AI/data/

GitHub push 为可选备份（需 token），非主要同步方式。
```

| 文件类型 | 本地 data/ | GitHub |
|----------|-----------|--------|
| `*.csv`, `*.smi`, `*.png`, `*.pdb` | ✅ Agent 自动写回 | 可选 push |
| `*.dcd`, `*.xtc`（MD 轨迹） | ❌ 太大 | 外部存储 |

## 三条流水线概览

| Notebook | 模态 | 核心工具链 | 文档 |
|----------|------|-----------|------|
| `OmniScreen_SM_Workflow` | 小分子 | RDKit, AutoDock Vina, OpenMM, MM/PBSA | [SM 模块说明](docs/workflows/SM_MODULES.md) |
| `OmniScreen_PE_Workflow` | 蛋白/多肽 | ESM-2, HDOCKlite, AlphaFold 3 | [PE 模块说明](docs/workflows/PE_MODULES.md) |
| `OmniScreen_NA_Workflow` | 核酸 | ViennaRNA, Bowtie2, AlphaFold 3 | [NA 模块说明](docs/workflows/NA_MODULES.md) |

## 数据说明

- 筛选结果 CSV、图表 PNG、小型 SMILES 库、受体 PDB → **自动 commit 到 GitHub**
- MD 轨迹（`.dcd`）、超大化合物库 → 不纳入 Git，见 `data/*/README.md`

## 方案局限与「漏检」

计算机辅助药物设计（AIDD）的核心目标不是「100% 不漏」，而是在数百万分子面前，用可承受的算力把候选范围尽快缩小。在这个过程中，**因软件计算出错、配置不当或算法近似而把好分子当成垃圾扔掉**，在工业界称为 **漏检（false negative）**。OmniScreen-AI 当前流水线同样存在这一权衡，主要体现在 Module 2–3 两个卡点。

### 卡点一：配体准备失败（`ligand_prep_failed`）

**原因：软件的「刻板印象」与能力边界。**

配体准备阶段，软件需要推测分子在生理环境（pH ≈ 7.4）下的质子化形态、电荷与可旋转键处理。部分新颖分子含有罕见杂环或复杂立体交叉结构时，经典拓扑工具（MGLTools、Open Babel 等）可能因参数库中缺乏先例，无法正确分配电荷或键参数，直接报错退出。

在本项目中，Module 3 使用 RDKit 3D 嵌入 + Open Babel 转 PDBQT；一旦准备失败，该分子不会进入 Vina 打分，在 `docking_scores.csv` 中标记为 `ligand_prep_failed`——**好药可能在此被无声淘汰**。

### 卡点二：Vina 对接失败或打分偏低（`vina_failed` / 低分淘汰）

**原因：刚性受体的「盲人摸象」与算力妥协。**

AutoDock Vina 为追求高通量速度，采用 **刚性受体对接**：把 PD-L1 当作固定不动的「锁」，只允许配体扭动构象去匹配。但真实蛋白是动态的：

- **诱导契合（Induced Fit）**：部分有效药物初接触时口袋尚小，静态对接会失败或给出极差分数（如 −4.0 kcal/mol）；若给予足够时间，配体可诱导口袋张开并形成稳定结合。Vina 的单次静态采样**看不到这一过程**。
- **搜索不充分**：默认 `exhaustiveness` 较低时，复杂口袋内可能未采样到合理姿态，导致运行失败或虚低分。

当前 OmniScreen SM 流程对通过初筛的分子按排序取 Top **250** 对接（`SM_CONFIG["max_dock"]`），进一步放大「排名靠后即被放弃」的漏检风险。

### 工业界常见应对策略

| 策略 | 思路 | 与本项目的关系 |
|------|------|----------------|
| **大基数放行** | Module 2 不过度收紧，让数千分子进入对接，用规模对抗规则漏检 | 当前库 ~2,300 条，初筛通过 ~864 条；可调宽 Lipinski 阈值或 `max_dock` |
| **深度学习对接并行** | DiffDock、Uni-Mol 等模型不依赖传统电荷参数，可捞回 Vina 失败分子 | 🔜 规划中，可作为 Module 3 补充分支 |
| **提高搜索精度** | 将 Vina `exhaustiveness` 从 8 提至 16–32，减少采样不足导致的失败 | 可在 Module 3 脚本中调整，代价是 CPU 时间成倍增加 |
| **动态复核** | MD / MM-PBSA 对幸存候选做诱导契合与结合稳定性验证 | Module 4–5（RunPod GPU）的设计初衷 |

### 设计哲学

> 宁可在大海里快速筛掉一大批，也不能让算力在第一步就瘫痪——但幸存下来的分子，仍需要一次「动态平反」。

因此 OmniScreen 采用 **Colab 漏斗（Module 0–3）+ RunPod 动力学（Module 4–5）** 的分层架构：前者负责速度与规模，后者负责对静态对接不确定性的补充验证。对接分数与初筛标签应视为 **排序线索**，而非最终成药性判决。

更细的模块级说明见 [SM 模块文档 — 局限性与假设](docs/workflows/SM_MODULES.md#15-当前局限性与假设)。

## License

MIT — 见 [LICENSE](LICENSE)
