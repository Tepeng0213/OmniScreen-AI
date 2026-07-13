# PD-L1 Receptor Structures

This directory holds receptor PDB files for screening (not tracked by Git).

## Recommended Structures

| PDB ID | Use | Notes |
|--------|------|------|
| `5N2F` / `6NM7` | Small-molecule line | PD-L1 dimer + small-molecule co-crystal; strip ligand for docking target |
| `4ZQK` | Protein/peptide line | PD-1/PD-L1 complex; defines binding interface hotspots |

## How to Obtain

```bash
# Example: download from RCSB
wget https://files.rcsb.org/download/5N2F.pdb -O data/receptor/5N2F.pdb
wget https://files.rcsb.org/download/4ZQK.pdb -O data/receptor/4ZQK.pdb
```

Or auto-download in each workflow notebook **Module 1**.
