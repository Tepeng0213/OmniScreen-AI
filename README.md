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
│   ├── OmniScreen_SM_Workflow.ipynb   # 小分子全流程
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
3. **配置 GitHub Token**（Colab 推送数据必需）：
   - 在 [GitHub Settings → Tokens](https://github.com/settings/tokens) 创建 Personal Access Token（勾选 `repo` 权限）
   - Colab 左侧 🔑 **Secrets** → 添加 `GITHUB_TOKEN`
4. 打开 `notebooks/OmniScreen_SM_Workflow.ipynb`，选择 **Colab** 内核
5. 按模块顺序运行 cell（Module 0 → Module 6）

> **算力建议**：Module 1–3 在 Colab CPU/GPU 免费层即可；Module 4–5（OpenMM / 自由能）建议在 RunPod 等 GPU 平台运行，详见各 notebook 内说明。

## 数据持久化策略

**所有运行数据写入 GitHub 仓库，绝不依赖 Colab 临时目录 `/content`。**

```
Colab 运行 Module 0
  → git clone https://github.com/Tepeng0213/OmniScreen-AI
  → 所有输出写入 /content/OmniScreen-AI/data/
  → 每个模块结束调用 persist_to_github()
  → git commit + push 到 GitHub

本地 Cursor
  → git pull 即可获取最新运行结果
```

| 文件类型 | 是否提交 GitHub |
|----------|----------------|
| `*.csv`, `*.smi`, `*.png` | ✅ 提交 |
| `*.pdb`（受体结构） | ✅ 提交 |
| `*.dcd`, `*.xtc`（MD 轨迹） | ❌ 太大，用外部存储 |

## 三条流水线概览

| Notebook | 模态 | 核心工具链 |
|----------|------|-----------|
| `OmniScreen_SM_Workflow` | 小分子 | RDKit, AutoDock Vina, OpenMM, MM/PBSA |
| `OmniScreen_PE_Workflow` | 蛋白/多肽 | ESM-2, HDOCKlite, AlphaFold 3 |
| `OmniScreen_NA_Workflow` | 核酸 | ViennaRNA, BWA-MEM2, AlphaFold 3 |

## 数据说明

- 筛选结果 CSV、图表 PNG、小型 SMILES 库、受体 PDB → **自动 commit 到 GitHub**
- MD 轨迹（`.dcd`）、超大化合物库 → 不纳入 Git，见 `data/*/README.md`

## License

MIT — 见 [LICENSE](LICENSE)
