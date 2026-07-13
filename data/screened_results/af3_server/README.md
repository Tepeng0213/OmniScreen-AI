# AlphaFold 3 Server Submission Materials

This directory is auto-organized from OmniScreen CPU-stage results (Modules 0–3) for upload to [AlphaFold Server](https://alphafoldserver.com).

> **Important**: The Server has **no** public API — you cannot drive it remotely from Cursor/Colab. Flow: notebook/this dir generates JSON → browser Upload JSON → download zip → back to notebook for visualization.

## Directory Layout

```
af3_server/
├── README.md                 # This guide
├── pe/                       # Protein/peptide: nanobody–PD-L1
│   ├── batch_top5.json       # ★ Recommended: upload 5 jobs at once
│   ├── PE_*.json             # Single-job backups
│   ├── manifest.csv          # Candidate metadata
│   └── pdl1_4zqk_chainB.fasta
├── na/                       # Nucleic acid: PD-L1 + siRNA
│   ├── batch_top5.json       # ★ Recommended: protein + antisense + sense (duplex)
│   ├── batch_top5_antisense_only.json  # Lighter: antisense only
│   ├── NA_*.json
│   └── manifest.csv
└── sm/                       # Small molecules: materials pack (Server capability limited)
    ├── SM_PDL1_5N2F_dimer_baseline.json  # PD-L1 dimer baseline (no custom ligand)
    ├── manifest_top5_ligands.csv
    ├── top5_ligands.smi      # For OpenMM MD / external docking
    └── pdl1_5n2f_chains.fasta
```

## Upload Steps (PE / NA)

1. Open https://alphafoldserver.com and sign in with Google  
2. Click **Upload JSON**  
3. Upload:
   - PE: `pe/batch_top5.json`
   - NA: `na/batch_top5.json` (if quota is tight, use `batch_top5_antisense_only.json`)
4. Check drafts in History → **Run** (uses daily free quota)  
5. When done, download zip and extract to e.g.:

```text
data/screened_results/af3_results/pe/
data/screened_results/af3_results/na/
```

Zip usually contains:
- `*_model.cif` / PDB
- `*_confidences.json` (ipTM / pTM / PAE)
- `*_job_request.json` (input for reproduction)

## Per-Line Notes

### PE (priority, best fit for Server)

| Item | Content |
|------|------|
| Complex | Top 5 nanobody mutants + PD-L1 (4ZQK chain B, 106 aa) |
| Source | Top 5 from `ppi_docking_scores.csv`: `CDR3_P98KV/A/L/I/P` |
| Expected | Antibody–antigen interface structure + ipTM / PAE |

### NA (good fit for Server)

| Item | Content |
|------|------|
| Complex | PD-L1 + siRNA antisense (+ sense duplex) |
| Source | Top 5 with `passed_offtarget` from `offtarget_filtered.csv` |
| Note | Biologically siRNA mainly targets **mRNA**; protein–nucleic jobs here are structural/spatial validation demos — interpret carefully |

### SM (Server barely supports custom ligands)

AlphaFold Server **only allows whitelist CCD ligands** (ATP, HEM, NAD, etc.), **not arbitrary SMILES**.  
Your top ligands (`CHEMBL19019_div`, etc.) cannot be submitted directly.

Prepared:

| File | Purpose |
|------|------|
| `SM_PDL1_5N2F_dimer_baseline.json` | Optional: run PD-L1 dimer baseline |
| `top5_ligands.smi` + `manifest_top5_ligands.csv` | For **OpenMM MD / MM-PBSA** (SM Modules 4–5 proper path) |
| `5N2F.pdb` (in repo) | MD docking receptor |

**SM recommended path**: skip custom ligand jobs on AF3 Server → use RunPod for OpenMM.

## Quota Suggestions

| Line | Jobs | Suggestion |
|------|--------|------|
| PE | 5 | Run first |
| NA | 5 (or antisense-only 5) | Next day or when quota allows |
| SM | 0–1 (dimer baseline only) | Optional |

## After Server Run — Back in Cursor

1. Extract zip to `data/screened_results/af3_results/{pe,na}/`  
2. In PE/NA notebook Module 4/6, parse CIF + confidence → 3D / radar plots  
3. `export_for_local_sync()` writes back locally (if visualization runs on Colab)

---

*Generated from OmniScreen CPU-stage results · 2026-07*
