"""Redraw paper Fig. F-G (peripheral actin, ancestor vs BasR G53V) from paper_stats.py output.

    python replication/plot_figure_FG.py stats_out/ [-o figure_FG.png]

Reads stats_out/well_medians_pooled.csv and stats_out/stats_pooled_lmm.csv
(written by paper_stats.py when given two or more experiments). Each point is
one well's median; the red line is the median of the uninfected-control wells;
stars are the BH-adjusted mixed-model p-values, as in the paper legend.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PANELS = [("cortical_intensity", "F", "Peripheral actin\nsignal intensity"),
          ("cortical_weighted_density", "G", "Peripheral actin\nbundle density")]
GROUPS = [("KpSen", "ancestor", "#9e9e9e"), ("BasR", "BasR G53V", "#f0a30a")]


def stars(q):
    for cut, s in ((1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (0.05, "*")):
        if q < cut:
            return s
    return "ns"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stats_dir")
    ap.add_argument("-o", "--out")
    args = ap.parse_args()
    d = Path(args.stats_dir)
    wm = pd.read_csv(d / "well_medians_pooled.csv")
    lmm = pd.read_csv(d / "stats_pooled_lmm.csv").set_index("feature")
    rng = np.random.default_rng(0)

    fig, axes = plt.subplots(2, 1, figsize=(3.2, 7.5))
    for ax, (col, letter, label) in zip(axes, PANELS):
        data = [wm.loc[wm.group == g, col].dropna().to_numpy() for g, _, _ in GROUPS]
        parts = ax.violinplot(data, showextrema=False, widths=0.8)
        for body, (_, _, color) in zip(parts["bodies"], GROUPS):
            body.set_facecolor(color)
            body.set_alpha(0.55)
        for i, (vals, (_, _, color)) in enumerate(zip(data, GROUPS), start=1):
            ax.boxplot(vals, positions=[i], widths=0.15, showfliers=False, medianprops={"color": "black"})
            ax.scatter(i + rng.uniform(-0.12, 0.12, len(vals)), vals, s=6, color=color, alpha=0.7, zorder=3)
        ctrl = wm.loc[wm.group == "Control", col].median()
        ax.axhline(ctrl, color="#d9534f", lw=1)
        top = max(v.max() for v in data)
        ax.plot([1, 1, 2, 2], [top * 1.04, top * 1.06, top * 1.06, top * 1.04], color="black", lw=0.8)
        ax.text(1.5, top * 1.07, stars(lmm.loc[col, "q_BH"]), ha="center")
        ax.set_xticks([1, 2], [lab for _, lab, _ in GROUPS], rotation=90)
        ax.set_ylabel(label)
        ax.set_title(letter, loc="left", fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    out = Path(args.out) if args.out else d / "figure_FG.png"
    fig.savefig(out, dpi=200)
    print("saved", out)


if __name__ == "__main__":
    main()
