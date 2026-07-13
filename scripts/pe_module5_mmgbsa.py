#!/usr/bin/env python
"""PE Module 5 — OpenMM MM-GBSA on AlphaFold 3 nanobody–target complexes.

Inputs
------
data/screened_results/af3_pe_metrics.csv
data/screened_results/af3_pe_complexes/<mut_id>_best.cif

Method
------
Amber14 + GBn2 implicit solvent (OpenMM CUDA):
  ΔG_bind ≈ E(complex) − E(receptor) − E(ligand)
Residue decomposition: pairwise VDW + ELE between nanobody residues and target.

Outputs
-------
data/screened_results/ppi_mmgbsa_summary.csv
data/screened_results/ppi_energy_decomposition.csv
data/screened_results/figures/fig_pe_mmgbsa_ranking.png
data/screened_results/figures/fig_pe_residue_energy_decomposition.png
data/screened_results/pe_complexes/<mut_id>_min.pdb
"""
from __future__ import annotations

import json
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1
from openmm import LangevinMiddleIntegrator, NonbondedForce, Platform
from openmm.app import (
    CutoffNonPeriodic,
    ForceField,
    HBonds,
    Modeller,
    PDBFile,
    Simulation,
)
from openmm.unit import (
    amu,
    elementary_charge,
    kelvin,
    kilojoule_per_mole,
    nanometer,
    picosecond,
    picoseconds,
)
from pdbfixer import PDBFixer
from tqdm import tqdm

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "screened_results"
COMPLEX_DIR = RESULTS / "pe_complexes"
FIG_DIR = RESULTS / "figures"
AF3_METRICS = RESULTS / "af3_pe_metrics.csv"
AF3_DIR = RESULTS / "af3_pe_complexes"

# AF3 PE jobs: chain A = nanobody, chain B = target
LIGAND_CHAIN = "A"
RECEPTOR_CHAIN = "B"

AA3_TO1 = {k.upper(): v.upper() for k, v in protein_letters_3to1.items()}


def prepare_system(cif_or_pdb: Path):
    fixer = PDBFixer(filename=str(cif_or_pdb))
    fixer.findMissingResidues()
    fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = ff.createSystem(
        fixer.topology,
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=2.0 * nanometer,
        constraints=HBonds,
        hydrogenMass=1.5 * amu,
    )
    return fixer, ff, system


def make_simulation(topology, system, positions):
    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    platform = Platform.getPlatformByName("CUDA")
    sim = Simulation(topology, system, integrator, platform, {"Precision": "mixed"})
    sim.context.setPositions(positions)
    return sim


def potential(sim) -> float:
    return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        kilojoule_per_mole
    )


def energy_subset(ff, topology, positions, keep_chains: set[str]) -> float:
    modeller = Modeller(topology, positions)
    delete = [r for r in modeller.topology.residues() if r.chain.id not in keep_chains]
    modeller.delete(delete)
    if modeller.topology.getNumAtoms() == 0:
        return 0.0
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=2.0 * nanometer,
        constraints=HBonds,
        hydrogenMass=1.5 * amu,
    )
    return potential(make_simulation(modeller.topology, system, modeller.positions))


def _get_nonbonded_force(system) -> NonbondedForce:
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, NonbondedForce):
            return force
    raise RuntimeError("No NonbondedForce in system")


def residue_pairwise_energy(system, topology, positions, chain_id: str, res_id: int, cutoff_nm: float = 2.0):
    other = RECEPTOR_CHAIN if chain_id == LIGAND_CHAIN else LIGAND_CHAIN
    nb = _get_nonbonded_force(system)
    xyz = np.array([p.value_in_unit(nanometer) for p in positions])
    params = []
    for i in range(nb.getNumParticles()):
        q, sigma, eps = nb.getParticleParameters(i)
        params.append(
            (
                q.value_in_unit(elementary_charge),
                sigma.value_in_unit(nanometer),
                eps.value_in_unit(kilojoule_per_mole),
            )
        )

    res_atoms, other_atoms = [], []
    for atom in topology.atoms():
        rid = int(atom.residue.id)
        cid = atom.residue.chain.id
        if cid == chain_id and rid == res_id:
            res_atoms.append(atom.index)
        elif cid == other:
            other_atoms.append(atom.index)
    if not res_atoms:
        return {"vdw_kJmol": np.nan, "elec_kJmol": np.nan, "total_kJmol": np.nan}

    k_e = 138.935456
    vdw = elec = 0.0
    cut2 = cutoff_nm * cutoff_nm
    for i in res_atoms:
        qi, si, ei = params[i]
        xi = xyz[i]
        for j in other_atoms:
            qj, sj, ej = params[j]
            d2 = float(np.sum((xi - xyz[j]) ** 2))
            if d2 > cut2 or d2 < 1e-12:
                continue
            r = np.sqrt(d2)
            sigma = 0.5 * (si + sj)
            epsilon = np.sqrt(max(ei * ej, 0.0))
            sr = sigma / r
            sr6 = sr ** 6
            vdw += 4.0 * epsilon * (sr6 * sr6 - sr6)
            elec += k_e * qi * qj / r
    return {"vdw_kJmol": vdw, "elec_kJmol": elec, "total_kJmol": vdw + elec}


def interface_residues(topology, positions, ligand_chain=LIGAND_CHAIN, receptor_chain=RECEPTOR_CHAIN, cutoff_nm=0.5):
    atom_pos = {a.index: np.array(positions[a.index].value_in_unit(nanometer)) for a in topology.atoms()}
    rec_xyz = []
    for res in topology.residues():
        if res.chain.id != receptor_chain:
            continue
        for atom in res.atoms():
            if atom.element.symbol != "H":
                rec_xyz.append(atom_pos[atom.index])
    if not rec_xyz:
        return []
    rec_xyz = np.asarray(rec_xyz)

    out = []
    for res in topology.residues():
        if res.chain.id != ligand_chain:
            continue
        lig = [atom_pos[a.index] for a in res.atoms() if a.element.symbol != "H"]
        if not lig:
            continue
        lig = np.asarray(lig)
        dmin = float(np.min(np.linalg.norm(lig[:, None, :] - rec_xyz[None, :, :], axis=2)))
        if dmin <= cutoff_nm:
            out.append((int(res.id), res.name, dmin))
    return out


def parse_mut_site(mut_id: str) -> int | None:
    # CDR3_P98KV -> 98
    import re

    m = re.match(r".*_P(\d+)[A-Z]{2}$", mut_id)
    return int(m.group(1)) if m else None


def analyze_complex(cif_path: Path, mut_id: str, meta: dict, minimize_steps: int = 200):
    fixer, ff, system = prepare_system(cif_path)
    chains = {c.id for c in fixer.topology.chains()}
    if LIGAND_CHAIN not in chains or RECEPTOR_CHAIN not in chains:
        raise RuntimeError(f"{mut_id}: expected chains A/B, found {chains}")

    sim = make_simulation(fixer.topology, system, fixer.positions)
    sim.minimizeEnergy(maxIterations=minimize_steps)
    state = sim.context.getState(getEnergy=True, getPositions=True)
    positions = state.getPositions()

    COMPLEX_DIR.mkdir(parents=True, exist_ok=True)
    min_pdb = COMPLEX_DIR / f"{mut_id}_min.pdb"
    with open(min_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, positions, f)

    e_c = potential(sim)
    e_r = energy_subset(ff, fixer.topology, positions, {RECEPTOR_CHAIN})
    e_l = energy_subset(ff, fixer.topology, positions, {LIGAND_CHAIN})
    dg = e_c - e_r - e_l

    mut_pos = parse_mut_site(mut_id)
    iface = interface_residues(fixer.topology, positions)
    if mut_pos is not None and not any(r[0] == mut_pos for r in iface):
        for res in fixer.topology.residues():
            if res.chain.id == LIGAND_CHAIN and int(res.id) == mut_pos:
                iface.append((mut_pos, res.name, np.nan))
                break

    rows = []
    for res_id, resname, dmin in iface:
        pair = residue_pairwise_energy(system, fixer.topology, positions, LIGAND_CHAIN, res_id)
        rows.append(
            {
                "mut_id": mut_id,
                "chain": LIGAND_CHAIN,
                "resnum": res_id,
                "resname": resname,
                "aa": AA3_TO1.get(resname, "X"),
                "is_mutated_site": bool(mut_pos is not None and res_id == mut_pos),
                "min_dist_nm": dmin,
                "E_vdw_kJmol": pair["vdw_kJmol"],
                "E_elec_kJmol": pair["elec_kJmol"],
                "E_residue_kJmol": pair["total_kJmol"],
            }
        )

    summary = {
        "mut_id": mut_id,
        "iptm": meta.get("iptm", np.nan),
        "ptm": meta.get("ptm", np.nan),
        "ranking_score": meta.get("ranking_score", np.nan),
        "esm_score": meta.get("esm_score", np.nan),
        "ppi_score": meta.get("ppi_score", np.nan),
        "E_complex_kJmol": e_c,
        "E_receptor_kJmol": e_r,
        "E_ligand_kJmol": e_l,
        "dG_bind_kJmol": dg,
        "dG_bind_kcalmol": dg / 4.184,
        "n_interface_residues": len(iface),
        "complex_cif": str(cif_path.relative_to(ROOT)),
        "minimized_pdb": str(min_pdb.relative_to(ROOT)),
    }
    return summary, pd.DataFrame(rows)


def plot_results(summary: pd.DataFrame, decomp: pd.DataFrame) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    s = summary.sort_values("dG_bind_kcalmol")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(range(len(s)), s["dG_bind_kcalmol"], color="#2a6f97")
    ax.set_xticks(range(len(s)))
    ax.set_xticklabels(s["mut_id"], rotation=40, ha="right")
    ax.set_ylabel("ΔG_bind (kcal/mol)")
    ax.set_title("PE Module 5 — MM-GBSA on AF3 complexes")
    ax.axhline(0, color="grey", lw=0.8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_pe_mmgbsa_ranking.png", dpi=180)
    plt.close(fig)

    top = s.iloc[0]["mut_id"]
    d = decomp[decomp["mut_id"] == top].sort_values("E_residue_kJmol")
    fig, ax = plt.subplots(figsize=(8, 5))
    labels = [f"{r.aa}{r.resnum}" for r in d.itertuples()]
    cols = ["#c1121f" if r.is_mutated_site else "#2a9d8f" for r in d.itertuples()]
    ax.barh(range(len(d)), d["E_residue_kJmol"] / 4.184, color=cols)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Residue VDW+ELE (kcal/mol)")
    ax.set_title(f"Residue energy decomposition — {top}")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_pe_residue_energy_decomposition.png", dpi=180)
    plt.close(fig)

    # Optional: ΔG vs ipTM
    if s["iptm"].notna().any():
        fig, ax = plt.subplots(figsize=(5.5, 4.5))
        ax.scatter(s["iptm"], s["dG_bind_kcalmol"], c="#264653", s=60)
        for r in s.itertuples():
            ax.annotate(r.mut_id.replace("CDR3_", ""), (r.iptm, r.dG_bind_kcalmol), fontsize=7, xytext=(4, 4), textcoords="offset points")
        ax.set_xlabel("AF3 ipTM")
        ax.set_ylabel("ΔG_bind (kcal/mol)")
        ax.set_title("MM-GBSA ΔG vs AF3 ipTM")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig_pe_mmgbsa_vs_iptm.png", dpi=180)
        plt.close(fig)


def main():
    if not AF3_METRICS.exists():
        raise FileNotFoundError(f"Missing {AF3_METRICS} — run PE Module 4 first")
    metrics = pd.read_csv(AF3_METRICS)
    if metrics.empty:
        raise RuntimeError("af3_pe_metrics.csv is empty")

    plats = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    print("OpenMM platforms:", plats)
    assert "CUDA" in plats, "CUDA required for PE Module 5"

    summaries, decomps = [], []
    print(f"Analyzing {len(metrics)} AF3 complexes...")
    for row in tqdm(metrics.itertuples(index=False), total=len(metrics), desc="MM-GBSA"):
        mut_id = row.mut_id
        cif = ROOT / row.complex_cif if not Path(row.complex_cif).is_absolute() else Path(row.complex_cif)
        if not cif.exists():
            alt = AF3_DIR / f"{mut_id}_best.cif"
            cif = alt
        if not cif.exists():
            print(f"[SKIP] missing CIF for {mut_id}")
            continue
        meta = row._asdict() if hasattr(row, "_asdict") else dict(zip(metrics.columns, row))
        try:
            summary, decomp = analyze_complex(cif, mut_id, meta)
        except Exception as exc:  # noqa: BLE001
            print(f"[FAIL] {mut_id}: {exc}")
            continue
        summaries.append(summary)
        decomps.append(decomp)
        print(f"  {mut_id}: ΔG = {summary['dG_bind_kcalmol']:.2f} kcal/mol | ipTM={meta.get('iptm')}")

    if not summaries:
        raise RuntimeError("All MM-GBSA jobs failed")

    summary_df = pd.DataFrame(summaries).sort_values("dG_bind_kcalmol")
    decomp_df = pd.concat(decomps, ignore_index=True)
    summary_path = RESULTS / "ppi_mmgbsa_summary.csv"
    decomp_path = RESULTS / "ppi_energy_decomposition.csv"
    summary_df.to_csv(summary_path, index=False)
    decomp_df.to_csv(decomp_path, index=False)
    plot_results(summary_df, decomp_df)

    meta_out = {
        "method": "OpenMM MM-GBSA (Amber14 + GBn2)",
        "inputs": "AF3 PE complexes (chain A nanobody / chain B target)",
        "platform": "CUDA",
        "n_complexes": int(len(summary_df)),
        "outputs": {
            "summary": str(summary_path.relative_to(ROOT)),
            "decomposition": str(decomp_path.relative_to(ROOT)),
        },
    }
    (RESULTS / "ppi_mmgbsa_meta.json").write_text(json.dumps(meta_out, indent=2))
    print("\nWrote", summary_path)
    print("Wrote", decomp_path)
    print(summary_df[["mut_id", "dG_bind_kcalmol", "iptm", "esm_score"]].to_string(index=False))


if __name__ == "__main__":
    os.chdir(ROOT)
    main()
