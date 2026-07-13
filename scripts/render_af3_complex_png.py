#!/usr/bin/env python
"""Render a static cartoon PNG from an AF3 complex CIF (PyMOL).

Used by PE/NA Module 4. Outputs structure only (no caption panel).
"""
from __future__ import annotations

import argparse
from pathlib import Path


def render_complex_png(
    cif_path: Path,
    out_png: Path,
    *,
    title: str = "",
    caption: str = "",
    chain_a_color: str = "skyblue",
    chain_b_color: str = "orange",
    chain_c_color: str | None = "forest",
    width: int = 1400,
    height: int = 1000,
) -> Path:
    """Render CIF to PNG. ``title`` / ``caption`` are kept for API compat but ignored."""
    del title, caption  # structure-only output
    from pymol import cmd

    cif_path = Path(cif_path)
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    cmd.reinitialize()
    cmd.load(str(cif_path), "cpx")
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.show("sticks", "cpx and (resn A+U+G+C+DA+DT+DG+DC+ADE+URA+GUA+CYT)")
    cmd.set("stick_radius", 0.15)
    cmd.color(chain_a_color, "cpx and chain A")
    cmd.color(chain_b_color, "cpx and chain B")
    if chain_c_color:
        cmd.color(chain_c_color, "cpx and chain C")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("antialias", 2)
    cmd.orient()
    cmd.zoom("cpx", 2)
    cmd.ray(width, height)
    cmd.png(str(out_png), dpi=150)
    return out_png


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="", help="Ignored (API compat)")
    ap.add_argument("--caption", default="", help="Ignored (API compat)")
    ap.add_argument("--no-chain-c", action="store_true", help="Skip coloring chain C")
    args = ap.parse_args()
    path = render_complex_png(
        Path(args.cif),
        Path(args.out),
        title=args.title,
        caption=args.caption,
        chain_c_color=None if args.no_chain_c else "forest",
    )
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
