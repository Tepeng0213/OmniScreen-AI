# Screening Results Output

Output files from each module of each workflow notebook live here.

**Persistence rule**: small results (CSV / PNG) are auto-committed via `persist_to_github()` to [GitHub](https://github.com/Tepeng0213/OmniScreen-AI). MD trajectories (`.dcd`) are too large for Git.

## Typical Small-Molecule Outputs

| File | Stage | Description |
|------|------|------|
| `chemical_space_props.csv` | Module 2 | Physicochemical properties + filter labels |
| `docking_scores.csv` | Module 3 | Vina docking score ranking |
| `top10_ligands.sdf` | Module 3 | Top 10 ligands (handoff to RunPod) |
| `md_rmsd.csv` | Module 4 | MD trajectory analysis |

## Typical Protein/Peptide Outputs

| File | Stage | Description |
|------|------|------|
| `mutation_scores.csv` | Module 2 | ESM-2 CDR mutation ΔLL scores |
| `ppi_docking_scores.csv` | Module 3 | Protein–protein interface scores |
| `pe_docking/*.pdb` | Module 3 | Nanobody structures (extended-CA or ESMFold) |
| `af3_pe_metrics.csv` | Module 4 | AF3 nanobody–PD-L1 ipTM / pTM |
| `af3_pe_complexes/*.cif` | Module 4 | AF3 best complex structures |
| `figures/fig5*.png`, `fig6*.png`, `fig_pe_af3_*` | Module 4/6 | CDR heatmap, PPI, AF3 ranking/3D |
| `af3_server/` | Module 4 prep | AlphaFold Server upload JSON (PE/NA/SM materials) |

See [PE modules](../../docs/workflows/PE_MODULES.md).

## Typical Nucleic-Acid Outputs

| File | Stage | Description |
|------|------|------|
| `CD274_mRNA.fasta` | Module 1 | PD-L1 mRNA sequence |
| `sirna_candidates.csv` | Module 2 | siRNA candidates + ViennaRNA scores |
| `offtarget_filtered.csv` | Module 3 | Candidates after off-target filtering |
| `sirna_offtarget.sam` | Module 3 | Bowtie2 alignment (chr22 demo) |
| `af3_na_metrics.csv` | Module 4 | AF3 protein–siRNA ipTM / pTM |
| `af3_na_complexes/*.cif` | Module 4 | AF3 best complex structures |
| `thermodynamics.csv` | Module 5 | Hybridization ΔG / hairpin ΔG |
| `figures/fig7*.png`, `fig8*.png`, `fig_na_*.png` | Module 4–6 | Efficacy, off-target, AF3, thermodynamics plots |

See [NA modules](../../docs/workflows/NA_MODULES.md).

## Cross-Platform Data Handoff

After Colab Modules 1–3, results are pushed to GitHub. On RunPod, `git clone` the same repo to continue Modules 4–5.

## Typical Protein/Peptide Outputs (continued)

| File | Stage | Description |
|------|------|------|
| `af3_pe_metrics.csv` / `af3_pe_complexes/` | Module 4 | AF3 complexes and ipTM |
| `ppi_mmgbsa_summary.csv` | Module 5 | MM-GBSA binding free energy ranking |
| `ppi_energy_decomposition.csv` | Module 5 | Per-residue VDW/ELE decomposition |
| `figures/fig_pe_mmgbsa_*.png` | Module 5 | ΔG ranking / residue contribution / vs ipTM |
