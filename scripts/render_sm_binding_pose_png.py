#!/usr/bin/env python
"""Render SM Top-1 protein–ligand binding-pose PNG via PyMOL."""
from __future__ import annotations

import argparse
from pathlib import Path


def render_binding_pose_png(
    complex_pdb: Path,
    out_png: Path,
    *,
    ligand_sele: str = "resn UNK",
    width: int = 1400,
    height: int = 1000,
) -> Path:
    from pymol import cmd

    complex_pdb = Path(complex_pdb)
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    cmd.reinitialize()
    cmd.load(str(complex_pdb), "cpx")
    cmd.hide("everything")
    cmd.remove("resn HOH+NA+CL")
    cmd.show("cartoon", "polymer")
    cmd.color("skyblue", "polymer and chain A")
    cmd.color("lightblue", "polymer and chain B")
    cmd.show("sticks", ligand_sele)
    cmd.color("orange", ligand_sele)
    cmd.set("stick_radius", 0.22)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("antialias", 2)
    cmd.orient(ligand_sele)
    cmd.zoom(ligand_sele, 12)
    cmd.ray(width, height)
    cmd.png(str(out_png), dpi=150)
    return out_png


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ligand-sele", default="resn UNK")
    args = ap.parse_args()
    path = render_binding_pose_png(Path(args.pdb), Path(args.out), ligand_sele=args.ligand_sele)
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
