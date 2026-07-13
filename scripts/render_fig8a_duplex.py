#!/usr/bin/env python
"""Render Fig 8a: Top-5 siRNA–mRNA duplex cartoons (structure-only visual)."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

NT_COLORS = {
    "A": "#e74c3c",
    "U": "#3498db",
    "G": "#f39c12",
    "C": "#27ae60",
    "T": "#3498db",
}


def draw_duplex_panel(ax, *, sirna_id: str, antisense: str, target: str, duplex_dg: float, efficacy: float, offtarget: int) -> None:
    """Draw antiparallel guide–target duplex as colored nucleotide beads + pair bars."""
    antisense = str(antisense).upper().replace("T", "U")
    target = str(target).upper().replace("T", "U")
    n = min(len(antisense), len(target))
    antisense = antisense[:n]
    # Display target 3'→5' under guide 5'→3'
    target_disp = target[:n][::-1]

    for i, nt in enumerate(antisense):
        ax.scatter(i, 1.15, s=220, c=NT_COLORS.get(nt, "#95a5a6"), zorder=3, edgecolors="#2c3e50", linewidths=0.4)
        ax.text(i, 1.15, nt, ha="center", va="center", fontsize=7.5, fontweight="bold", color="white", zorder=4)
    for i, nt in enumerate(target_disp):
        ax.scatter(i, 0.0, s=220, c=NT_COLORS.get(nt, "#95a5a6"), zorder=3, edgecolors="#2c3e50", linewidths=0.4)
        ax.text(i, 0.0, nt, ha="center", va="center", fontsize=7.5, fontweight="bold", color="white", zorder=4)

    for i in range(n):
        a, b = antisense[i], target_disp[i]
        paired = {a, b} in ({"A", "U"}, {"G", "C"}, {"G", "U"})
        ax.plot([i, i], [0.22, 0.93], color="#2c3e50" if paired else "#bdc3c7", lw=2.0 if paired else 1.0, ls="-" if paired else ":", zorder=1)

    ax.text(-1.1, 1.15, "5'", fontsize=8, ha="center", va="center")
    ax.text(n + 0.1, 1.15, "3'", fontsize=8, ha="center", va="center")
    ax.text(-1.1, 0.0, "3'", fontsize=8, ha="center", va="center")
    ax.text(n + 0.1, 0.0, "5'", fontsize=8, ha="center", va="center")
    ax.text(-1.8, 1.15, "guide", fontsize=7, ha="right", va="center", color="#34495e")
    ax.text(-1.8, 0.0, "target", fontsize=7, ha="right", va="center", color="#34495e")

    ax.set_xlim(-2.4, n + 0.8)
    ax.set_ylim(-0.45, 1.75)
    ax.axis("off")
    ax.set_title(
        f"{sirna_id}   efficacy={efficacy:.1f}   off-target={int(offtarget)}   duplex ΔG={duplex_dg:.1f} kcal/mol",
        fontsize=9,
        loc="left",
        pad=2,
    )


def render_fig8a(
    *,
    offtarget_csv: Path,
    thermo_csv: Path,
    out_png: Path,
    top_n: int = 5,
) -> Path:
    df_off = pd.read_csv(offtarget_csv).head(top_n)
    df_th = pd.read_csv(thermo_csv).set_index("sirna_id")

    fig, axes = plt.subplots(len(df_off), 1, figsize=(12.5, 1.85 * len(df_off)))
    if len(df_off) == 1:
        axes = [axes]

    for ax, (_, row) in zip(axes, df_off.iterrows()):
        sid = row["sirna_id"]
        th = df_th.loc[sid]
        draw_duplex_panel(
            ax,
            sirna_id=sid,
            antisense=th["antisense_seq"],
            target=th["target_seq"],
            duplex_dg=float(th["duplex_dg"]),
            efficacy=float(row["efficacy_score"]),
            offtarget=int(row["offtarget_hits_chr22"]),
        )

    # nucleotide legend on last panel
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markeredgecolor="#2c3e50", markersize=8, label=nt)
        for nt, c in [("A", NT_COLORS["A"]), ("U", NT_COLORS["U"]), ("G", NT_COLORS["G"]), ("C", NT_COLORS["C"])]
    ]
    axes[-1].legend(handles=handles, loc="lower right", ncol=4, frameon=False, fontsize=8, bbox_to_anchor=(1.0, -0.15))

    fig.suptitle("Fig 8a — Top-5 siRNA–mRNA Duplex (guide 5′→3′ / target 3′→5′)", y=0.995, fontsize=12)
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=[0, 0.02, 1, 0.98])
    fig.savefig(out_png, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_png


def main():
    root = Path(__file__).resolve().parents[1]
    out = render_fig8a(
        offtarget_csv=root / "data/screened_results/offtarget_filtered.csv",
        thermo_csv=root / "data/screened_results/thermodynamics.csv",
        out_png=root / "data/screened_results/figures/fig8a_rna_structure.png",
    )
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
