#!/usr/bin/env python
"""SM Module 5 — MM/GBSA free-energy analysis from Module 4 MD trajectories."""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1
from openff.toolkit import Molecule
from openff.units import unit as offunit
from openmm import LangevinMiddleIntegrator, NonbondedForce, Platform, Vec3
from openmm.app import ForceField, HBonds, Modeller, NoCutoff, Simulation
from openmm.unit import amu, elementary_charge, kelvin, kilojoule_per_mole, nanometer, picosecond, picoseconds
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from pdbfixer import PDBFixer
from tqdm import tqdm

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "screened_results"
MD_ROOT = RESULTS / "md"
FIG_DIR = RESULTS / "figures"
AA3_TO1 = {k.upper(): v.upper() for k, v in protein_letters_3to1.items()}


def make_sim(topology, system, positions):
    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    platform = Platform.getPlatformByName("CUDA")
    sim = Simulation(topology, system, integrator, platform, {"Precision": "mixed"})
    sim.context.setPositions(positions)
    return sim


def potential(sim) -> float:
    return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole)


def prepare_receptor():
    fixer = PDBFixer(filename=str(ROOT / "data/receptor/5N2F.pdb"))
    fixer.findMissingResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    return fixer


def _get_nb(system) -> NonbondedForce:
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, NonbondedForce):
            return force
    raise RuntimeError("No NonbondedForce")


def residue_pair_energies(system, topology, positions, n_protein_atoms: int, cutoff_nm: float = 1.2):
    nb = _get_nb(system)
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

    lig_atoms = list(range(n_protein_atoms, topology.getNumAtoms()))
    res_atoms: dict[tuple, list[int]] = {}
    for atom in topology.atoms():
        if atom.index >= n_protein_atoms:
            continue
        key = (atom.residue.chain.id, int(atom.residue.id), atom.residue.name)
        res_atoms.setdefault(key, []).append(atom.index)

    k_e = 138.935456
    cut2 = cutoff_nm * cutoff_nm
    rows = []
    for (chain, resnum, resname), atoms in res_atoms.items():
        vdw = elec = 0.0
        dmin = 1e9
        for i in atoms:
            qi, si, ei = params[i]
            xi = xyz[i]
            for j in lig_atoms:
                d2 = float(np.sum((xi - xyz[j]) ** 2))
                dmin = min(dmin, d2)
                if d2 > cut2 or d2 < 1e-12:
                    continue
                r = np.sqrt(d2)
                qj, sj, ej = params[j]
                sigma = 0.5 * (si + sj)
                epsilon = np.sqrt(max(ei * ej, 0.0))
                sr = sigma / r
                sr6 = sr ** 6
                vdw += 4.0 * epsilon * (sr6 * sr6 - sr6)
                elec += k_e * qi * qj / r
        if dmin > cut2:
            continue
        rows.append(
            {
                "chain": chain,
                "resnum": resnum,
                "resname": resname,
                "aa": AA3_TO1.get(resname.upper(), "X"),
                "min_dist_nm": float(np.sqrt(dmin)),
                "E_vdw_kJmol": vdw,
                "E_elec_kJmol": elec,
                "E_residue_kJmol": vdw + elec,
            }
        )
    return pd.DataFrame(rows)


def analyze_ligand(mol_id: str, receptor: PDBFixer, n_frames: int):
    out_dir = MD_ROOT / mol_id
    meta = json.loads((out_dir / "md_meta.json").read_text())
    traj = md.load(str(out_dir / "trajectory.dcd"), top=str(out_dir / "topology.pdb"))

    n_prot = int(meta["n_protein_atoms"])
    n_lig = int(meta["n_ligand_atoms"])
    prot_slice = slice(0, n_prot)
    lig_slice = slice(n_prot, n_prot + n_lig)

    offmol = Molecule.from_smiles(meta["smiles"], allow_undefined_stereo=True)
    offmol.generate_conformers(n_conformers=1)
    try:
        offmol.assign_partial_charges("gasteiger")
    except Exception:
        offmol.assign_partial_charges("mmff94")
    if offmol.n_atoms != n_lig:
        raise RuntimeError(f"{mol_id}: OpenFF atoms {offmol.n_atoms} != MD ligand atoms {n_lig}")

    smirnoff = SMIRNOFFTemplateGenerator(molecules=offmol, forcefield="openff-2.1.0")
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    ff.registerTemplateGenerator(smirnoff.generator)

    idxs = np.linspace(0, len(traj) - 1, num=min(n_frames, len(traj)), dtype=int)
    frame_rows = []
    decomp_df = None

    for idx in idxs:
        xyz = traj.xyz[idx]  # nm
        prot_xyz = xyz[prot_slice]
        lig_xyz = xyz[lig_slice]
        if len(prot_xyz) != receptor.topology.getNumAtoms():
            raise RuntimeError(
                f"{mol_id}: protein atoms {len(prot_xyz)} != receptor {receptor.topology.getNumAtoms()}"
            )

        offmol.conformers[0] = lig_xyz * offunit.nanometer
        lig_top = offmol.to_topology().to_openmm()
        lig_pos = offmol.conformers[0].to_openmm()
        prot_pos = [Vec3(*map(float, row)) for row in prot_xyz] * nanometer

        complex_mod = Modeller(receptor.topology, prot_pos)
        complex_mod.add(lig_top, lig_pos)

        system = ff.createSystem(
            complex_mod.topology,
            nonbondedMethod=NoCutoff,
            constraints=HBonds,
            hydrogenMass=1.5 * amu,
        )
        e_c = potential(make_sim(complex_mod.topology, system, complex_mod.positions))

        sys_r = ff.createSystem(
            receptor.topology,
            nonbondedMethod=NoCutoff,
            constraints=HBonds,
            hydrogenMass=1.5 * amu,
        )
        e_r = potential(make_sim(receptor.topology, sys_r, prot_pos))

        sys_l = ff.createSystem(
            lig_top,
            nonbondedMethod=NoCutoff,
            constraints=HBonds,
            hydrogenMass=1.5 * amu,
        )
        e_l = potential(make_sim(lig_top, sys_l, lig_pos))
        dg = e_c - e_r - e_l
        frame_rows.append(
            {
                "mol_id": mol_id,
                "frame": int(idx),
                "E_complex_kJmol": e_c,
                "E_receptor_kJmol": e_r,
                "E_ligand_kJmol": e_l,
                "dG_bind_kJmol": dg,
                "dG_bind_kcalmol": dg / 4.184,
            }
        )
        if decomp_df is None:
            decomp = residue_pair_energies(system, complex_mod.topology, complex_mod.positions, n_prot)
            decomp.insert(0, "mol_id", mol_id)
            decomp_df = decomp

    frame_df = pd.DataFrame(frame_rows)
    summary = {
        "mol_id": mol_id,
        "smiles": meta["smiles"],
        "vina_score": meta["vina_score"],
        "ligand_rmsd_mean_A": meta.get("ligand_rmsd_mean_A", np.nan),
        "n_frames": len(frame_df),
        "dG_bind_kcalmol_mean": float(frame_df["dG_bind_kcalmol"].mean()),
        "dG_bind_kcalmol_std": float(frame_df["dG_bind_kcalmol"].std(ddof=0)),
        "dG_bind_kJmol_mean": float(frame_df["dG_bind_kJmol"].mean()),
    }
    return summary, frame_df, decomp_df if decomp_df is not None else pd.DataFrame()


def plot_results(summary: pd.DataFrame, decomp: pd.DataFrame) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    s = summary.sort_values("dG_bind_kcalmol_mean")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(range(len(s)), s["dG_bind_kcalmol_mean"], yerr=s["dG_bind_kcalmol_std"], color="#1d3557", capsize=3)
    ax.set_xticks(range(len(s)))
    ax.set_xticklabels(s["mol_id"], rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("ΔG_bind (kcal/mol)")
    ax.set_title("SM Module 5 — MM/GBSA from MD frames")
    ax.axhline(0, color="grey", lw=0.8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_sm_mmpbsa_ranking.png", dpi=180)
    plt.close(fig)

    top = s.iloc[0]["mol_id"]
    d = decomp[decomp["mol_id"] == top].nsmallest(20, "E_residue_kJmol")
    fig, ax = plt.subplots(figsize=(8, 5))
    labels = [f"{r.aa}{r.resnum}{r.chain}" for r in d.itertuples()]
    ax.barh(range(len(d)), d["E_residue_kJmol"] / 4.184, color="#457b9d")
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Residue VDW+ELE (kcal/mol)")
    ax.set_title(f"Residue contributions — {top}")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_sm_mmpbsa_residue.png", dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-frames", type=int, default=5)
    args = ap.parse_args()

    plats = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    print("OpenMM platforms:", plats)
    assert "CUDA" in plats

    md_summary = RESULTS / "md_rmsd.csv"
    if not md_summary.exists():
        raise FileNotFoundError("Run Module 4 first (md_rmsd.csv missing)")
    ligands = pd.read_csv(md_summary)["mol_id"].tolist()
    receptor = prepare_receptor()

    summaries, frames, decomps = [], [], []
    for mol_id in tqdm(ligands, desc="MM/GBSA"):
        try:
            summary, frame_df, decomp_df = analyze_ligand(mol_id, receptor, args.n_frames)
        except Exception as exc:  # noqa: BLE001
            print(f"[FAIL] {mol_id}: {exc}")
            continue
        summaries.append(summary)
        frames.append(frame_df)
        decomps.append(decomp_df)
        print(
            f"  {mol_id}: ΔG = {summary['dG_bind_kcalmol_mean']:.2f} ± "
            f"{summary['dG_bind_kcalmol_std']:.2f} kcal/mol"
        )

    if not summaries:
        raise RuntimeError("All MM/GBSA jobs failed")

    summary_df = pd.DataFrame(summaries).sort_values("dG_bind_kcalmol_mean")
    frame_df = pd.concat(frames, ignore_index=True)
    decomp_df = pd.concat(decomps, ignore_index=True)
    summary_df.to_csv(RESULTS / "mmpbsa_results.csv", index=False)
    frame_df.to_csv(RESULTS / "mmpbsa_frame_energies.csv", index=False)
    decomp_df.to_csv(RESULTS / "mmpbsa_residue_decomposition.csv", index=False)
    plot_results(summary_df, decomp_df)
    print("\nWrote", RESULTS / "mmpbsa_results.csv")
    print(
        summary_df[
            ["mol_id", "vina_score", "dG_bind_kcalmol_mean", "dG_bind_kcalmol_std", "ligand_rmsd_mean_A"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
