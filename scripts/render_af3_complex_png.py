#!/usr/bin/env python
"""Render a static cartoon PNG from an AF3 complex CIF (PyMOL).

Used by PE/NA Module 4 to replace text-only caption cards.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def render_complex_png(
    cif_path: Path,
    out_png: Path,
    *,
    title: str,
    caption: str,
    chain_a_color: str = "skyblue",
    chain_b_color: str = "orange",
    width: int = 1400,
    height: int = 1000,
) -> Path:
    from pymol import cmd
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import image as mpimg

    cif_path = Path(cif_path)
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    raw = out_png.with_suffix(".raw.png")

    cmd.reinitialize()
    cmd.load(str(cif_path), "cpx")
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.color(chain_a_color, "cpx and chain A")
    cmd.color(chain_b_color, "cpx and chain B")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("antialias", 2)
    cmd.orient()
    cmd.zoom("cpx", 2)
    cmd.ray(width, height)
    cmd.png(str(raw), dpi=150)

    img = mpimg.imread(raw)
    fig = plt.figure(figsize=(9.5, 7.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[5.8, 1.15], hspace=0.08)
    ax = fig.add_subplot(gs[0])
    ax.imshow(img)
    ax.axis("off")
    ax.set_title(title, fontsize=13, pad=8)

    ax2 = fig.add_subplot(gs[1])
    ax2.axis("off")
    ax2.text(0.01, 0.5, caption, fontsize=10, family="DejaVu Sans", va="center")
    fig.savefig(out_png, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    try:
        raw.unlink()
    except OSError:
        pass
    return out_png


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--caption", required=True)
    args = ap.parse_args()
    path = render_complex_png(Path(args.cif), Path(args.out), title=args.title, caption=args.caption)
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
