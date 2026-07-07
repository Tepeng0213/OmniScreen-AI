# PD-L1 受体结构

本目录存放筛选用的受体 PDB 文件（不纳入 Git 追踪）。

## 推荐结构

| PDB ID | 用途 | 说明 |
|--------|------|------|
| `5N2F` / `6NM7` | 小分子线 | PD-L1 二聚体 + 小分子共晶，剥离配体后作对接靶点 |
| `4ZQK` | 蛋白/多肽线 | PD-1/PD-L1 复合物，定义结合界面热点 |

## 获取方式

```bash
# 示例：从 RCSB 下载
wget https://files.rcsb.org/download/5N2F.pdb -O data/receptor/5N2F.pdb
wget https://files.rcsb.org/download/4ZQK.pdb -O data/receptor/4ZQK.pdb
```

或在各 workflow notebook 的 **Module 1** 中自动下载。
