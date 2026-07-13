#!/usr/bin/env python
"""Regenerate all PE/SM/NA figures with English titles (no MD/docking recompute)."""
from __future__ import annotations

import json
import os
import sys
import traceback
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_nb(name: str) -> dict:
    return json.loads((ROOT / "notebooks" / name).read_text(encoding="utf-8"))


def cell_src(cell: dict) -> str:
    return "".join(cell.get("source", []))


def strip_magics(src: str) -> str:
    lines = []
    for line in src.splitlines():
        s = line.strip()
        if s.startswith("%") or s.startswith("!"):
            continue
        lines.append(line)
    return "\n".join(lines)


def exec_src(src: str, ns: dict, label: str) -> bool:
    code = strip_magics(src)
    if not code.strip():
        return True
    try:
        exec(compile(code, label, "exec"), ns, ns)
        print(f"  OK  {label}")
        return True
    except Exception as e:
        print(f"  FAIL {label}: {type(e).__name__}: {e}")
        traceback.print_exc(limit=2)
        return False


def base_ns() -> dict:
    ns: dict = {"__name__": "__main__"}
    exec(
        f"""
from pathlib import Path
import os, sys, json, warnings, subprocess, shutil, math, re, time
warnings.filterwarnings("ignore")
PROJECT_ROOT = Path(r"{ROOT}")
os.chdir(PROJECT_ROOT)
sys.path.insert(0, str(PROJECT_ROOT))
PATHS = {{
    "project_root": PROJECT_ROOT,
    "receptor": PROJECT_ROOT / "data" / "receptor",
    "raw": PROJECT_ROOT / "data" / "raw_libraries",
    "results": PROJECT_ROOT / "data" / "screened_results",
}}
def export_for_local_sync(files):
    print(f"(skip sync) {{len(list(files))}} files")
def persist_to_github(*a, **k):
    print("(skip github persist)")
""",
        ns,
        ns,
    )
    return ns


def regenerate_sm() -> None:
    print("=== SM ===")
    nb = load_nb("OmniScreen_SM_Workflow.ipynb")
    ns = base_ns()
    # Module 0 for SM_CONFIG + helpers (skip clone-heavy parts by providing stubs if needed)
    m0 = cell_src(nb["cells"][4])
    exec_src(m0, ns, "SM:Module0")
    ns["RECEPTOR_PDB"] = "5N2F"
    # Module 6.0–6.7
    for i in (21, 23, 25, 27, 29, 31, 33, 35):
        exec_src(cell_src(nb["cells"][i]), ns, f"SM:cell{i}")

    FIG = ROOT / "data/screened_results/figures"
    # RMSD curves from per-ligand csv
    rows = []
    for p in sorted((ROOT / "data/screened_results/md").glob("*/rmsd.csv")):
        rows.append(pd.read_csv(p))
    if rows:
        df = pd.concat(rows, ignore_index=True)
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for mol, g in df.groupby("mol_id"):
            ax.plot(g["time_ps"], g["ligand_rmsd_A"], label=mol)
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel("Ligand RMSD (Å)")
        ax.set_title("SM Module 4 — Ligand RMSD (demo MD)")
        ax.legend(fontsize=7, loc="best")
        fig.tight_layout()
        fig.savefig(FIG / "fig_sm_md_rmsd.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
        print("  OK  fig_sm_md_rmsd.png")

    mm = ROOT / "data/screened_results/mmpbsa_results.csv"
    if mm.exists():
        df = pd.read_csv(mm).sort_values("dG_bind_kcalmol_mean")
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.barh(df["mol_id"], df["dG_bind_kcalmol_mean"], color="#3498db")
        ax.set_xlabel("ΔG_bind (kcal/mol)")
        ax.set_title("SM Module 5 — MM/GBSA from MD frames")
        fig.tight_layout()
        fig.savefig(FIG / "fig_sm_mmpbsa_ranking.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
        print("  OK  fig_sm_mmpbsa_ranking.png")

    res = ROOT / "data/screened_results/mmpbsa_residue_decomposition.csv"
    if res.exists():
        df = pd.read_csv(res)
        # pick best (most negative) mol from mmpbsa_results if available
        top = df["mol_id"].iloc[0]
        if mm.exists():
            top = pd.read_csv(mm).sort_values("dG_bind_kcalmol_mean").iloc[0]["mol_id"]
        sub = df[df["mol_id"] == top].copy()
        sub["label"] = sub["resname"].astype(str) + sub["resnum"].astype(str)
        sub = sub.nsmallest(20, "E_residue_kJmol")
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.barh(sub["label"], sub["E_residue_kJmol"], color="#e67e22")
        ax.set_xlabel("E_residue (kJ/mol)")
        ax.set_title(f"Residue contributions — {top}")
        fig.tight_layout()
        fig.savefig(FIG / "fig_sm_mmpbsa_residue.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
        print("  OK  fig_sm_mmpbsa_residue.png")

    # Binding pose PyMOL
    try:
        from importlib.util import module_from_spec, spec_from_file_location

        spec = spec_from_file_location("rbp", ROOT / "scripts/render_sm_binding_pose_png.py")
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        pdb = ROOT / "data/screened_results/md/CHEMBL19019_div/complex_min.pdb"
        if pdb.exists():
            mod.render_binding_pose_png(pdb, FIG / "fig_3d_binding_pose.png")
            print("  OK  fig_3d_binding_pose.png (PyMOL)")
    except Exception as e:
        print(f"  FAIL binding pose: {e}")


def regenerate_pe() -> None:
    print("=== PE ===")
    nb = load_nb("OmniScreen_PE_Workflow.ipynb")
    ns = base_ns()
    # Minimal Module 0 if needed
    for i, c in enumerate(nb["cells"]):
        src = cell_src(c)
        if c.get("cell_type") == "code" and src.strip().startswith("# @title Module 0"):
            exec_src(src, ns, "PE:Module0")
            break
    # Module 6
    exec_src(cell_src(nb["cells"][14]), ns, "PE:Module6")

    FIG = ROOT / "data/screened_results/figures"
    # AF3 ipTM ranking
    metrics = ROOT / "data/screened_results/af3_pe_metrics.csv"
    if metrics.exists():
        df = pd.read_csv(metrics).sort_values("iptm", ascending=False)
        fig, ax = plt.subplots(figsize=(8, 4.5))
        colors = ["#e67e22" if i == 0 else "#3498db" for i in range(len(df))]
        ax.barh(df["mut_id"][::-1], df["iptm"][::-1], color=colors[::-1], edgecolor="white")
        ax.axvline(0.6, color="#e74c3c", ls="--", lw=1, label="ipTM≈0.6 guideline")
        ax.set_xlabel("ipTM")
        ax.set_title("Fig PE-AF3 — Nanobody–PD-L1 interface confidence")
        ax.legend(loc="lower right")
        fig.tight_layout()
        fig.savefig(FIG / "fig_pe_af3_iptm_ranking.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("  OK  fig_pe_af3_iptm_ranking.png")

        top = df.iloc[0]["mut_id"]
        cif = ROOT / f"data/screened_results/af3_pe_complexes/{top}_best.cif"
        if cif.exists():
            from importlib.util import module_from_spec, spec_from_file_location

            spec = spec_from_file_location("raf3", ROOT / "scripts/render_af3_complex_png.py")
            mod = module_from_spec(spec)
            spec.loader.exec_module(mod)
            mod.render_complex_png(
                cif,
                FIG / "fig_pe_af3_complex.png",
                title="",
                caption="",
                chain_c_color=None,
            )
            print("  OK  fig_pe_af3_complex.png")

    # MM-GBSA plots from CSV (English titles)
    summ = ROOT / "data/screened_results/ppi_mmgbsa_summary.csv"
    if summ.exists():
        df = pd.read_csv(summ).sort_values("dG_bind_kcalmol")
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.barh(df["mut_id"], df["dG_bind_kcalmol"], color="#3498db")
        ax.set_xlabel("ΔG_bind (kcal/mol)")
        ax.set_title("PE Module 5 — MM-GBSA on AF3 complexes")
        fig.tight_layout()
        fig.savefig(FIG / "fig_pe_mmgbsa_ranking.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
        print("  OK  fig_pe_mmgbsa_ranking.png")

        if "iptm" in df.columns or metrics.exists():
            m = pd.read_csv(metrics) if metrics.exists() else df
            merged = df.merge(m[["mut_id", "iptm"]], on="mut_id", how="left", suffixes=("", "_m"))
            if "iptm" not in merged.columns:
                merged["iptm"] = merged.get("iptm_m")
            fig, ax = plt.subplots(figsize=(6.5, 5))
            ax.scatter(merged["iptm"], merged["dG_bind_kcalmol"], s=80, c="#e67e22")
            for _, r in merged.iterrows():
                ax.annotate(r["mut_id"], (r["iptm"], r["dG_bind_kcalmol"]), fontsize=7, xytext=(4, 4), textcoords="offset points")
            ax.set_xlabel("AF3 ipTM")
            ax.set_ylabel("ΔG_bind (kcal/mol)")
            ax.set_title("MM-GBSA ΔG vs AF3 ipTM")
            fig.tight_layout()
            fig.savefig(FIG / "fig_pe_mmgbsa_vs_iptm.png", dpi=180, bbox_inches="tight")
            plt.close(fig)
            print("  OK  fig_pe_mmgbsa_vs_iptm.png")

    decomp = ROOT / "data/screened_results/ppi_energy_decomposition.csv"
    if decomp.exists() and summ.exists():
        top = pd.read_csv(summ).sort_values("dG_bind_kcalmol").iloc[0]["mut_id"]
        df = pd.read_csv(decomp)
        sub = df[df["mut_id"] == top].copy() if "mut_id" in df.columns else df.copy()
        # flexible columns
        if "resname" in sub.columns and "resnum" in sub.columns:
            sub["label"] = sub["resname"].astype(str) + sub["resnum"].astype(str)
        else:
            sub["label"] = sub.iloc[:, 1].astype(str)
        ecol = "E_total_kcalmol" if "E_total_kcalmol" in sub.columns else [c for c in sub.columns if "kcal" in c.lower() or c.startswith("E_")][-1]
        sub = sub.nsmallest(20, ecol)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.barh(sub["label"], sub[ecol], color="#e67e22")
        ax.set_xlabel(ecol)
        ax.set_title(f"Residue energy decomposition — {top}")
        fig.tight_layout()
        fig.savefig(FIG / "fig_pe_residue_energy_decomposition.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
        print("  OK  fig_pe_residue_energy_decomposition.png")


def regenerate_na() -> None:
    print("=== NA ===")
    nb = load_nb("OmniScreen_NA_Workflow.ipynb")
    ns = base_ns()
    for i, c in enumerate(nb["cells"]):
        src = cell_src(c)
        if c.get("cell_type") == "code" and src.strip().startswith("# @title Module 0"):
            exec_src(src, ns, "NA:Module0")
            break
    exec_src(cell_src(nb["cells"][16]), ns, "NA:Module6")

    FIG = ROOT / "data/screened_results/figures"
    # AF3 ranking + complex
    metrics = ROOT / "data/screened_results/af3_na_metrics.csv"
    if metrics.exists():
        df = pd.read_csv(metrics).sort_values("iptm", ascending=False)
        fig, ax = plt.subplots(figsize=(8, 4.5))
        colors = ["#e67e22" if i == 0 else "#3498db" for i in range(len(df))]
        ax.barh(df["sirna_id"][::-1], df["iptm"][::-1], color=colors[::-1], edgecolor="white")
        ax.axvline(0.6, color="#e74c3c", ls="--", lw=1, label="ipTM≈0.6 guideline")
        ax.set_xlabel("ipTM")
        ax.set_title("Fig NA-AF3 — Protein–siRNA complex interface confidence")
        ax.legend(loc="lower right")
        fig.tight_layout()
        fig.savefig(FIG / "fig_na_af3_iptm_ranking.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("  OK  fig_na_af3_iptm_ranking.png")

        top = df.iloc[0]["sirna_id"]
        cif = ROOT / f"data/screened_results/af3_na_complexes/{top}_best.cif"
        if cif.exists():
            from importlib.util import module_from_spec, spec_from_file_location

            spec = spec_from_file_location("raf3", ROOT / "scripts/render_af3_complex_png.py")
            mod = module_from_spec(spec)
            spec.loader.exec_module(mod)
            mod.render_complex_png(cif, FIG / "fig_na_af3_complex.png", title="", caption="")
            print("  OK  fig_na_af3_complex.png")

    # thermo scatter
    thermo = ROOT / "data/screened_results/thermodynamics.csv"
    if thermo.exists():
        df = pd.read_csv(thermo)
        fig, ax = plt.subplots(figsize=(7, 5))
        passed = df["passed_thermo"] if "passed_thermo" in df.columns else pd.Series([True] * len(df))
        ax.scatter(df.loc[~passed, "duplex_dg"], df.loc[~passed, "hairpin_dg"], c="#bdc3c7", s=20, label="fail", alpha=0.6)
        ax.scatter(df.loc[passed, "duplex_dg"], df.loc[passed, "hairpin_dg"], c="#2ecc71", s=20, label="pass", alpha=0.7)
        n_pass = int(passed.sum())
        ax.set_xlabel("Duplex ΔG (kcal/mol)")
        ax.set_ylabel("Hairpin ΔG (kcal/mol)")
        ax.set_title(f"Fig NA-Thermo — Hybridization vs self-fold ({n_pass}/{len(df)} pass)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(FIG / "fig_na_thermo_scatter.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("  OK  fig_na_thermo_scatter.png")

    # fig8a duplex
    try:
        from importlib.util import module_from_spec, spec_from_file_location

        spec = spec_from_file_location("f8a", ROOT / "scripts/render_fig8a_duplex.py")
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        mod.render_fig8a(
            offtarget_csv=ROOT / "data/screened_results/offtarget_filtered.csv",
            thermo_csv=ROOT / "data/screened_results/thermodynamics.csv",
            out_png=FIG / "fig8a_rna_structure.png",
        )
        print("  OK  fig8a_rna_structure.png")
    except Exception as e:
        print(f"  FAIL fig8a: {e}")


def main():
    regenerate_sm()
    regenerate_pe()
    regenerate_na()
    print("=== DONE ===")


if __name__ == "__main__":
    main()
