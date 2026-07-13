# OmniScreen-AI

**OmniScreen-AI: A Multi-Modal High-Throughput Drug Screening & Dynamic Verification Platform for PD-L1**

A full-modality (small molecules / macrocyclic peptides / nucleic acids) high-throughput screening and dynamic kinetics verification platform for the tumor immunotherapy target **PD-L1**.

## Project Highlights

- **Three-modality parallel funnel**: chemical space (small molecules), sequence space (antibodies/peptides), RNA space (siRNA/Aptamer)
- **Four-stage closed-loop validation**: AI pre-screening → structural docking → MD kinetics → free energy calculation
- **Notebook-driven**: all logic lives in 3 workflow notebooks, modular blocks, one-click Colab reproduction
- **Tiered compute**: Colab CPU (pre-screen + docking) → RunPod GPU (MD + free energy)

## Repository Structure

```
OmniScreen-AI/
├── notebooks/
│   ├── OmniScreen_SM_Workflow.ipynb          # Small-molecule full pipeline (ChEMBL real compound library)
│   ├── OmniScreen_PE_Workflow.ipynb   # Protein/peptide full pipeline
│   └── OmniScreen_NA_Workflow.ipynb   # Nucleic-acid drug full pipeline
├── data/
│   ├── receptor/          # PD-L1 receptor PDB (5N2F, 4ZQK, etc.)
│   ├── raw_libraries/     # Initial candidate libraries (SMILES / FASTA / mRNA)
│   └── screened_results/  # Stage-wise screening outputs
├── docker/                # Optional: standardized runtime environment
├── pyproject.toml
└── requirements.txt
```

## Quick Start (Cursor + Colab)

1. Clone this repo and open it in Cursor
2. Install the [Google Colab](https://marketplace.visualstudio.com/items?itemName=Google.colab) and [Notebook MCP](https://marketplace.visualstudio.com/items?itemName=olavovieiradecarvalho.notebook-mcp-server) extensions
3. **Configure GitHub Token** (required for pushing data):
   - Create a PAT at [GitHub Settings → Tokens](https://github.com/settings/tokens) (check `repo` scope)
   - **Cursor + Colab extension**: run the "Configure GitHub Token" cell in the notebook and paste the token manually (Colab Secrets only work in the web UI; the Cursor kernel cannot read them)
   - **Colab web UI**: add `GITHUB_TOKEN` in 🔑 Secrets; code will read it automatically
4. Open `notebooks/OmniScreen_SM_Workflow.ipynb` and select the **Colab** kernel
5. Run cells in module order (Module 0 → Module 6)

> **Compute tips**: Modules 1–3 run fine on Colab free-tier CPU/GPU; Modules 4–5 (OpenMM / free energy) are best on GPU platforms like RunPod — see in-notebook notes.

## Data Persistence Strategy (Colab + Local Auto-Sync)

**Colab runs in the cloud and cannot write directly to local disk.** We use "cloud compute → pack outputs → Agent writes back locally":

```
Colab runs modules
  → export_for_local_sync() base64-encodes data/ files into cell output
  → Cursor Agent parses sync markers
  → Writes automatically to local /Users/schmit/Documents/OmniScreen-AI/data/

GitHub push is optional backup (needs token), not the primary sync path.
```

| File type | Local data/ | GitHub |
|----------|-----------|--------|
| `*.csv`, `*.smi`, `*.png`, `*.pdb` | ✅ Agent auto-write | Optional push |
| `*.dcd`, `*.xtc` (MD trajectories) | ❌ Too large | External storage |

## Three Pipeline Overview

| Notebook | Modality | Core toolchain | Docs |
|----------|------|-----------|------|
| `OmniScreen_SM_Workflow` | Small molecules | RDKit, AutoDock Vina, OpenMM, MM/PBSA | [SM modules](docs/workflows/SM_MODULES.md) |
| `OmniScreen_PE_Workflow` | Protein/peptide | ESM-2, HDOCKlite, AlphaFold 3 | [PE modules](docs/workflows/PE_MODULES.md) |
| `OmniScreen_NA_Workflow` | Nucleic acids | ViennaRNA, Bowtie2, AlphaFold 3 | [NA modules](docs/workflows/NA_MODULES.md) |

## Data Notes

- Screening CSVs, chart PNGs, small SMILES libraries, receptor PDB → **auto-commit to GitHub**
- MD trajectories (`.dcd`), oversized compound libraries → not in Git; see `data/*/README.md`

## Limitations & False Negatives

The core goal of computer-aided drug design (AIDD) is not "100% recall" but to shrink the candidate space with affordable compute among millions of molecules. **Discarding good molecules due to software errors, misconfiguration, or algorithmic approximations** is called a **false negative** in industry. The OmniScreen-AI pipeline has the same trade-off, mainly at Modules 2–3.

### Bottleneck 1: Ligand preparation failure (`ligand_prep_failed`)

**Cause: rigid software assumptions and capability limits.**

During ligand prep, tools infer protonation state, charge, and rotatable bonds at physiological pH (~7.4). Novel molecules with rare heterocycles or complex stereochemistry may lack precedents in classical topology tools (MGLTools, Open Babel, etc.), causing charge or parameter assignment to fail.

In this project, Module 3 uses RDKit 3D embedding + Open Babel → PDBQT; on failure the molecule never enters Vina scoring and is marked `ligand_prep_failed` in `docking_scores.csv` — **promising drugs may be silently dropped here**.

### Bottleneck 2: Vina docking failure or low scores (`vina_failed` / low-score elimination)

**Cause: rigid receptor "blind probing" and compute compromises.**

AutoDock Vina uses **rigid receptor docking** for throughput: PD-L1 is a fixed "lock" and only the ligand rotates. Real proteins are dynamic:

- **Induced fit**: some actives bind when the pocket is still small; static docking fails or scores poorly (e.g. −4.0 kcal/mol); given time, the ligand can open the pocket and stabilize. Vina's single static sample **misses this**.
- **Insufficient search**: low default `exhaustiveness` may miss reasonable poses in complex pockets, causing failure or artificially low scores.

The current OmniScreen SM flow docks only the top **250** after pre-screen (`SM_CONFIG["max_dock"]`), increasing false-negative risk for lower-ranked molecules.

### Common industry mitigations

| Strategy | Idea | Relation to this project |
|------|------|----------------|
| **Large funnel throughput** | Don't over-tighten Module 2; let thousands dock; use scale against rule-based misses | Current library ~2,300; ~864 pass pre-screen; widen Lipinski or `max_dock` |
| **Deep-learning docking in parallel** | DiffDock, Uni-Mol, etc. avoid classical charge params; recover Vina failures | 🔜 Planned as Module 3 supplement |
| **Higher search precision** | Raise Vina `exhaustiveness` from 8 to 16–32 | Adjustable in Module 3; CPU time grows |
| **Dynamic re-validation** | MD / MM-PBSA for induced fit and binding stability | Modules 4–5 (RunPod GPU) purpose |

### Design philosophy

> Better to quickly filter a huge pool than stall compute at step one — but survivors still need dynamic "appeal."

OmniScreen uses **Colab funnel (Modules 0–3) + RunPod kinetics (Modules 4–5)**: speed and scale first, then validation of static docking uncertainty. Docking scores and pre-screen labels are **ranking hints**, not final druggability verdicts.

See [SM modules — limitations & assumptions](docs/workflows/SM_MODULES.md#15-current-limitations-and-assumptions) for module-level detail.

## License

MIT — see [LICENSE](LICENSE)
