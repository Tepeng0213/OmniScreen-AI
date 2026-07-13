#!/usr/bin/env python
"""SM Module 4 — OpenMM explicit-solvent MD for top Vina ligands.

Inputs
------
data/screened_results/docking_scores.csv
data/receptor/5N2F.pdb   (pocket defined by co-crystal ligand 8HW)

Outputs
-------
data/screened_results/top10_ligands.sdf
data/screened_results/md/<mol_id>/complex_min.pdb
data/screened_results/md/<mol_id>/trajectory.dcd
data/screened_results/md/<mol_id>/topology.pdb
data/screened_results/md_rmsd.csv
data/screened_results/figures/fig_sm_md_rmsd.png
"""
from __future__ import annotations

import argparse
import json
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from openff.toolkit import Molecule
from openff.units import unit as offunit
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform
from openmm.app import (
    HBonds,
    Modeller,
    PDBFile,
    PME,
    Simulation,
    StateDataReporter,
    DCDReporter,
)
from openmm.unit import amu, atmosphere, kelvin, molar, nanometer, picosecond, picoseconds
from openmmforcefields.generators import SystemGenerator
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "screened_results"
RECEPTOR = ROOT / "data" / "receptor" / "5N2F.pdb"
MD_ROOT = RESULTS / "md"
FIG_DIR = RESULTS / "figures"


def pocket_com_nm(pdb_path: Path) -> np.ndarray:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rec", str(pdb_path))
    coords = []
    for res in structure[0]["B"]:
        if res.resname == "8HW":
            for atom in res:
                if atom.element != "H":
                    coords.append(atom.coord)
    if not coords:
        # fallback: protein COM
        for atom in structure.get_atoms():
            if atom.element != "H":
                coords.append(atom.coord)
    return np.asarray(coords).mean(axis=0) / 10.0  # Å → nm


def write_top10_sdf(df: pd.DataFrame, out_sdf: Path, n: int = 10) -> pd.DataFrame:
    top = df[df["status"] == "ok"].sort_values("vina_score").head(n).copy()
    writer = Chem.SDWriter(str(out_sdf))
    keep_rows = []
    for _, row in top.iterrows():
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) != 0:
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass
        mol.SetProp("_Name", str(row["mol_id"]))
        mol.SetProp("vina_score", f"{row['vina_score']:.3f}")
        writer.write(mol)
        keep_rows.append(row)
    writer.close()
    return pd.DataFrame(keep_rows)


def build_ligand(smiles: str, pocket: np.ndarray) -> Molecule:
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol is None:
        raise ValueError(f"Bad SMILES: {smiles}")
    rdmol = Chem.AddHs(rdmol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(rdmol, params) != 0:
        if AllChem.EmbedMolecule(rdmol, AllChem.ETKDG()) != 0:
            raise RuntimeError("3D embed failed")
    try:
        AllChem.MMFFOptimizeMolecule(rdmol, maxIters=200)
    except Exception:
        pass
    offmol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
    xyz = np.array(offmol.conformers[0].m_as("nanometer"))
    xyz = xyz - xyz.mean(axis=0) + pocket
    offmol.conformers[0] = xyz * offunit.nanometer
    # Gasteiger is far faster than AM1-BCC/sqm for demo throughput
    try:
        offmol.assign_partial_charges("gasteiger")
    except Exception:
        offmol.assign_partial_charges("mmff94")
    return offmol


def prepare_receptor(pdb_path: Path) -> PDBFixer:
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    return fixer


def run_md_for_ligand(
    mol_id: str,
    smiles: str,
    vina_score: float,
    pocket: np.ndarray,
    receptor: PDBFixer,
    nvt_ps: float,
    npt_ps: float,
    prod_ps: float,
    report_interval_ps: float,
) -> dict:
    out_dir = MD_ROOT / mol_id
    out_dir.mkdir(parents=True, exist_ok=True)

    offmol = build_ligand(smiles, pocket)
    ligand_top = offmol.to_topology().to_openmm()
    ligand_pos = offmol.conformers[0].to_openmm()

    modeller = Modeller(receptor.topology, receptor.positions)
    n_protein_atoms = modeller.topology.getNumAtoms()
    modeller.add(ligand_top, ligand_pos)
    n_complex_atoms = modeller.topology.getNumAtoms()
    n_ligand_atoms = n_complex_atoms - n_protein_atoms

    system_generator = SystemGenerator(
        forcefields=["amber14-all.xml", "amber14/tip3pfb.xml"],
        small_molecule_forcefield="openff-2.1.0",
        molecules=[offmol],
        forcefield_kwargs={
            "constraints": HBonds,
            "rigidWater": True,
            "removeCMMotion": True,
            "hydrogenMass": 1.5 * amu,
        },
        periodic_forcefield_kwargs={"nonbondedMethod": PME, "nonbondedCutoff": 1.0 * nanometer},
    )
    modeller.addSolvent(system_generator.forcefield, padding=1.0 * nanometer, ionicStrength=0.15 * molar)
    system = system_generator.create_system(modeller.topology)
    system.addForce(MonteCarloBarostat(1 * atmosphere, 300 * kelvin))

    dt = 0.002  # ps
    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, dt * picoseconds)
    platform = Platform.getPlatformByName("CUDA")
    sim = Simulation(modeller.topology, system, integrator, platform, {"Precision": "mixed"})
    sim.context.setPositions(modeller.positions)

    print(f"  [{mol_id}] atoms={modeller.topology.getNumAtoms()} minimize...")
    sim.minimizeEnergy(maxIterations=200)
    with open(out_dir / "complex_min.pdb", "w") as f:
        PDBFile.writeFile(
            sim.topology,
            sim.context.getState(getPositions=True).getPositions(),
            f,
        )
    # dry topology for later analysis (no water/ions) — save full top for DCD
    with open(out_dir / "topology.pdb", "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)

    sim.context.setVelocitiesToTemperature(300 * kelvin)
    nvt_steps = int(nvt_ps / dt)
    npt_steps = int(npt_ps / dt)
    prod_steps = int(prod_ps / dt)
    report_steps = max(1, int(report_interval_ps / dt))

    print(f"  [{mol_id}] NVT {nvt_ps} ps...")
    sim.step(nvt_steps)
    print(f"  [{mol_id}] NPT {npt_ps} ps...")
    sim.step(npt_steps)

    dcd_path = out_dir / "trajectory.dcd"
    log_path = out_dir / "md.log"
    # clear old reporters
    sim.reporters.clear()
    sim.reporters.append(DCDReporter(str(dcd_path), report_steps))
    sim.reporters.append(
        StateDataReporter(
            str(log_path),
            report_steps,
            step=True,
            time=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
            speed=True,
        )
    )
    print(f"  [{mol_id}] production {prod_ps} ps...")
    sim.step(prod_steps)

    # RMSD: ligand heavy atoms vs first production frame (atom index range)
    traj = md.load(str(dcd_path), top=str(out_dir / "topology.pdb"))
    lig_atoms = []
    for atom in traj.topology.atoms:
        if n_protein_atoms <= atom.index < n_protein_atoms + n_ligand_atoms and atom.element.symbol != "H":
            lig_atoms.append(atom.index)
    lig_atoms = np.array(lig_atoms, dtype=int)
    if len(lig_atoms) == 0:
        lig_atoms = np.arange(n_protein_atoms, n_protein_atoms + n_ligand_atoms)

    rmsd = md.rmsd(traj, traj, 0, atom_indices=lig_atoms) * 10.0  # nm→Å
    times = np.arange(len(rmsd)) * report_interval_ps
    rmsd_df = pd.DataFrame({"time_ps": times, "ligand_rmsd_A": rmsd, "mol_id": mol_id})
    rmsd_df.to_csv(out_dir / "rmsd.csv", index=False)

    meta = {
        "mol_id": mol_id,
        "smiles": smiles,
        "vina_score": float(vina_score),
        "n_atoms_solvated": int(modeller.topology.getNumAtoms()),
        "n_protein_atoms": int(n_protein_atoms),
        "n_ligand_atoms": int(n_ligand_atoms),
        "protein_atom_start": 0,
        "ligand_atom_start": int(n_protein_atoms),
        "nvt_ps": nvt_ps,
        "npt_ps": npt_ps,
        "prod_ps": prod_ps,
        "n_frames": int(len(rmsd)),
        "ligand_rmsd_mean_A": float(np.mean(rmsd)),
        "ligand_rmsd_final_A": float(rmsd[-1]),
        "ligand_rmsd_max_A": float(np.max(rmsd)),
        "trajectory": str((out_dir / "trajectory.dcd").relative_to(ROOT)),
        "topology": str((out_dir / "topology.pdb").relative_to(ROOT)),
    }
    (out_dir / "md_meta.json").write_text(json.dumps(meta, indent=2))
    return meta, rmsd_df


def plot_rmsd(all_rmsd: pd.DataFrame) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for mol_id, g in all_rmsd.groupby("mol_id"):
        ax.plot(g["time_ps"], g["ligand_rmsd_A"], label=mol_id, lw=1.5)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Ligand RMSD (Å)")
    ax.set_title("SM Module 4 — Ligand RMSD (demo MD)")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_sm_md_rmsd.png", dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-n", type=int, default=5, help="Number of top Vina ligands")
    ap.add_argument("--nvt-ps", type=float, default=25.0)
    ap.add_argument("--npt-ps", type=float, default=25.0)
    ap.add_argument("--prod-ps", type=float, default=200.0)
    ap.add_argument("--report-ps", type=float, default=10.0)
    args = ap.parse_args()

    plats = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    print("OpenMM platforms:", plats)
    assert "CUDA" in plats

    dock = pd.read_csv(RESULTS / "docking_scores.csv")
    top = write_top10_sdf(dock, RESULTS / "top10_ligands.sdf", n=max(10, args.top_n))
    run_df = top.head(args.top_n)
    print(f"Top ligands for MD:\n{run_df[['mol_id','vina_score']].to_string(index=False)}")

    pocket = pocket_com_nm(RECEPTOR)
    print("Pocket COM (nm):", pocket)
    receptor = prepare_receptor(RECEPTOR)

    metas, rmsds = [], []
    for row in tqdm(run_df.itertuples(index=False), total=len(run_df), desc="MD"):
        done_meta = MD_ROOT / row.mol_id / "md_meta.json"
        if done_meta.exists() and not os.environ.get("FORCE_RERUN"):
            print(f"  [{row.mol_id}] skip existing")
            meta = json.loads(done_meta.read_text())
            rmsd_df = pd.read_csv(MD_ROOT / row.mol_id / "rmsd.csv")
            metas.append(meta)
            rmsds.append(rmsd_df)
            continue
        try:
            meta, rmsd_df = run_md_for_ligand(
                row.mol_id,
                row.smiles,
                float(row.vina_score),
                pocket,
                receptor,
                args.nvt_ps,
                args.npt_ps,
                args.prod_ps,
                args.report_ps,
            )
        except Exception as exc:  # noqa: BLE001
            print(f"[FAIL] {row.mol_id}: {exc}")
            continue
        metas.append(meta)
        rmsds.append(rmsd_df)

    if not metas:
        raise RuntimeError("All MD jobs failed")

    summary = pd.DataFrame(metas)
    all_rmsd = pd.concat(rmsds, ignore_index=True)
    summary_path = RESULTS / "md_rmsd.csv"
    # store per-ligand summary + keep detailed trajectories in md/
    summary.to_csv(summary_path, index=False)
    all_rmsd.to_csv(RESULTS / "md_rmsd_timeseries.csv", index=False)
    plot_rmsd(all_rmsd)
    print("\nWrote", summary_path)
    print(summary[["mol_id", "vina_score", "ligand_rmsd_mean_A", "ligand_rmsd_final_A"]].to_string(index=False))


if __name__ == "__main__":
    main()
