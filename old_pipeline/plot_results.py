"""
plot_results.py — Visualize TopoActin V3 pipeline CSV output.

Usage (no arguments): processes all EXPERIMENT_CONFIGS defined below.
    python plot_results.py

Usage (single file): backwards-compatible, one-off single-file mode.
    python plot_results.py <path_to_csv> [output_name]

Each experiment produces a multi-page PDF containing:
  - Plate map (optional, toggled by SHOW_PLATE_MAP)
  - Summary heatmap(s) of z-scored group means
  - One page per metric category
      single mode: violin + strip plots
      paired mode: grouped boxplots (one box per condition per strain)
  - Ranked metric changes chart  (paired mode)
  - PCA overview                 (requires scikit-learn)
  - Statistical significance      (paired mode; requires scipy)
  - Metric correlation heatmap(s)
"""

import sys
import pathlib
import warnings
warnings.filterwarnings("ignore", message=".*tight_layout.*")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.backends.backend_pdf as pdf_backend
import seaborn as sns


# ---------------------------------------------------------------------------
# Metric category definitions  (edit column lists to customise output)
# ---------------------------------------------------------------------------

METRIC_GROUPS = {
    "Cell Morphology": [
        "cell_area",
        "cell_solidity",
        "cell_eccentricity",
        "cell_perimeter",
    ],
    "Skeleton Structure": [
        "total_length",
        "num_branches",
        "mean_branch_length",
        "num_junctions",
        "num_endpoints",
        "skeleton_density",
        "mean_tortuosity_edge",
    ],
    "Intensity & Weighting": [
        "total_network_intensity",
        "mean_bundle_intensity",
        "weighted_total_length",
        # intensity_capture_ratio removed — identical to fidelity
        "fidelity",
    ],
    "Zonal Actin Density": [
        "inner_weighted_density",
        "mid_weighted_density",
        "cortical_weighted_density",
        "cortical_inner_ratio_weighted",
        "cortical_mid_ratio_weighted",
        "mid_inner_ratio_weighted",
    ],
    "Radial Actin Profile": [
        "zone_1_density",
        "zone_2_density",
        "zone_3_density",
        # inner_density / mid_density / outer_density removed — identical to zone_N_density
        "peripheral_gradient",
    ],
    "Advanced Topology": [
        "zone_1_lacunarity",
        "zone_3_lacunarity",
        "spectral_entropy",
    ],
    "Bacterial Burden": [
        "intracellular_bact_area",
        "area_ratio",
        "bact_density_perinuclear",
        "bact_density_middle",
        "bact_density_cortical",
    ],
}


# ---------------------------------------------------------------------------
# Experiment configuration  (edit this block to switch experiments)
# ---------------------------------------------------------------------------
#
# mode "single":
#   {"mode": "single", "csv": "path/to/file.csv", "output_name": "optional_name"}
#
# mode "paired":
#   {"mode": "paired",
#    "csv1": "path/condition_A.csv", "label1": "1h",
#    "csv2": "path/condition_B.csv", "label2": "2h",
#    "output_name": "A549_mScarlet_1h_vs_2h"}
#
# Run all: python plot_results.py
# One-off: python plot_results.py <file.csv> [output_name]
# ---------------------------------------------------------------------------

_BASE = pathlib.Path(__file__).parent   # directory containing this script

EXPERIMENT_CONFIGS = [
    {
        "mode": "paired",
        "csv1": str(_BASE / "cytoskeleton_analysis_v3_A549_mScarlet_1h_2026-02-10_20260227_053200.csv"),
        "label1": "1h",
        "csv2": str(_BASE / "cytoskeleton_analysis_v3_A549_mScarlet_2h_2026-02-10_20260226_235412.csv"),
        "label2": "2h",
        "output_name": "A549_mScarlet_1h_vs_2h",
    },
    {
        "mode": "paired",
        "csv1": str(_BASE / "cytoskeleton_analysis_v3_A549_mScarlet_2026-02-12_20260226_152258.csv"),
        "label1": "Control",
        "csv2": str(_BASE / "cytoskeleton_analysis_v3_A549_cytD_mScarlet_2026-02-12_20260226_210738.csv"),
        "label2": "+CytD 2.5\u00b5L",
        "output_name": "A549_mScarlet_Control_vs_CytD_2.5uL",
    },
    {
        # Combined view: both experiments side by side per metric, shared y-axis per row
        # so effect sizes are directly comparable between the two experiments.
        "mode": "combined",
        "output_name": "A549_combined_timepoint_vs_cytD",
        "pairs": [
            {
                "csv1": str(_BASE / "cytoskeleton_analysis_v3_A549_mScarlet_1h_2026-02-10_20260227_053200.csv"),
                "label1": "1h",
                "csv2": str(_BASE / "cytoskeleton_analysis_v3_A549_mScarlet_2h_2026-02-10_20260226_235412.csv"),
                "label2": "2h",
                "comparison_label": "Timepoint  (1h vs 2h)",
            },
            {
                "csv1": str(_BASE / "cytoskeleton_analysis_v3_A549_mScarlet_2026-02-12_20260226_152258.csv"),
                "label1": "Control",
                "csv2": str(_BASE / "cytoskeleton_analysis_v3_A549_cytD_mScarlet_2026-02-12_20260226_210738.csv"),
                "label2": "+CytD 2.5\u00b5L",
                "comparison_label": "CytD Treatment  (Control vs +CytD 2.5\u00b5L)",
            },
        ],
    },
]


# ---------------------------------------------------------------------------
# Plate Map Configuration  (edit to match your physical plate layout)
# ---------------------------------------------------------------------------
#
# SHOW_PLATE_MAP = True  → adds a plate layout diagram as the first PDF page
# SHOW_PLATE_MAP = False → skips the plate map entirely
#
# PLATE_LAYOUT is a 2D list that mirrors the physical plate grid:
#   outer list = rows  (A, B, C, ...)  — one entry per row
#   inner list = cols  (1, 2, 3, ...)  — one entry per column
#
#   Cell values:
#     "GroupName"   → well assigned to that treatment group
#     "+GroupName"  → same group but CytD-treated (rendered with a hatch pattern)
#     None          → empty well
#
# Resize for other formats:
#   24-well: PLATE_ROWS = list("ABCD"),  PLATE_COLS = list(range(1, 7))
#    6-well: PLATE_ROWS = list("AB"),    PLATE_COLS = list(range(1, 4))
# ---------------------------------------------------------------------------

SHOW_PLATE_MAP = True

PLATE_ROWS = list("ABCDEFGH")    # 8 rows  (96-well plate)
PLATE_COLS = list(range(1, 13))  # 12 cols (96-well plate)

PLATE_LAYOUT = [
    # Col:    1           2           3         4         5          6          7      8      9        10       11      12
    ["Control",  "Control",  "KpSen",  "KpSen",  "Nissle",  "Nissle",  "BW",  "BW",  "LF82",  "LF82",  "C+M",  "Cpom"],  # Row A
    ["+Control", "+Control", "+KpSen", "+KpSen", "+Nissle", "+Nissle", "+BW", "+BW", "+LF82", "+LF82", "+C+M", "+Cpom"],  # Row B (CytD)
    [None, None, None, None, None, None, None, None, None, None, None, None],  # Row C
    [None, None, None, None, None, None, None, None, None, None, None, None],  # Row D
    [None, None, None, None, None, None, None, None, None, None, None, None],  # Row E
    [None, None, None, None, None, None, None, None, None, None, None, None],  # Row F
    [None, None, None, None, None, None, None, None, None, None, None, None],  # Row G
    [None, None, None, None, None, None, None, None, None, None, None, None],  # Row H
]


# ---------------------------------------------------------------------------
# Global settings  (tweak these without touching the plotting functions)
# ---------------------------------------------------------------------------

# Preferred group display order — unlisted groups are appended alphabetically
GROUP_ORDER_PREF = ["Control", "KpSen", "Nissle", "BW", "LF82", "C+M", "Cpom"]

# Group colour palette (tab10 gives 10 distinct colours)
PALETTE = sns.color_palette("tab10")

# Two-condition palette for paired plots (muted is more colorblind-friendly)
_COND_COLORS = sns.color_palette("muted")

# IQR multiplier for outlier removal (standard = 1.5; increase to be less aggressive)
IQR_K = 1.5


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _group_order(df: pd.DataFrame) -> list:
    """Return groups in preferred order; unknown groups appended alphabetically."""
    present = set(df["group"].dropna().unique())
    ordered = [g for g in GROUP_ORDER_PREF if g in present]
    extra = sorted(present - set(ordered))
    return ordered + extra


def _make_palette(groups: list) -> dict:
    return {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(groups)}


def _nice_label(col: str) -> str:
    return col.replace("_", " ").title()


def _all_metric_cols(df: pd.DataFrame) -> list:
    """Flatten METRIC_GROUPS to a deduplicated, ordered list of columns present in df."""
    seen, out = set(), []
    for cols in METRIC_GROUPS.values():
        for c in cols:
            if c not in seen and c in df.columns:
                seen.add(c)
                out.append(c)
    return out


def _drop_constant_cols(df: pd.DataFrame, cols: list, min_unique: int = 2) -> list:
    """Return only cols with >= min_unique distinct non-NaN values; print dropped names."""
    kept, dropped = [], []
    for c in cols:
        if c not in df.columns:
            continue
        if df[c].nunique(dropna=True) >= min_unique:
            kept.append(c)
        else:
            dropped.append(c)
    if dropped:
        print(f"  [auto-drop] constant/single-value columns removed from plots: {dropped}")
    return kept


def _remove_outliers_iqr(df: pd.DataFrame, cols: list, k: float = IQR_K) -> pd.DataFrame:
    """Replace values outside [Q1 - k*IQR, Q3 + k*IQR] per column per group with NaN."""
    df = df.copy()
    for col in cols:
        if col not in df.columns or not pd.api.types.is_numeric_dtype(df[col]):
            continue
        for grp_name in df["group"].dropna().unique():
            mask_grp = df["group"] == grp_name
            vals = df.loc[mask_grp, col]
            q1, q3 = vals.quantile(0.25), vals.quantile(0.75)
            iqr = q3 - q1
            lo, hi = q1 - k * iqr, q3 + k * iqr
            outlier = mask_grp & ((df[col] < lo) | (df[col] > hi))
            df.loc[outlier, col] = np.nan
    return df


# ---------------------------------------------------------------------------
# Plate map page
# ---------------------------------------------------------------------------

def plot_plate_map_page(groups: list, palette: dict) -> plt.Figure:
    """Render a plate layout diagram; CytD-treated wells are hatched."""
    nrows = len(PLATE_ROWS)
    ncols = len(PLATE_COLS)
    fig_w = ncols * 1.1 + 3.5
    fig_h = nrows * 1.1 + 2.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for ri in range(nrows):
        for ci in range(ncols):
            if ri < len(PLATE_LAYOUT) and ci < len(PLATE_LAYOUT[ri]):
                well_val = PLATE_LAYOUT[ri][ci]
            else:
                well_val = None

            if well_val is None:
                color = "#e0e0e0"
                hatch = None
                text = ""
            elif well_val.startswith("+"):
                base_grp = well_val[1:]
                color = palette.get(base_grp, "#cccccc")
                hatch = "///"
                text = base_grp
            else:
                color = palette.get(well_val, "#cccccc")
                hatch = None
                text = well_val

            rect = mpatches.FancyBboxPatch(
                (ci + 0.05, nrows - ri - 0.95),
                0.88, 0.88,
                boxstyle="round,pad=0.04",
                facecolor=color,
                hatch=hatch,
                linewidth=0.8,
                edgecolor="white",
            )
            ax.add_patch(rect)

            if text:
                ax.text(
                    ci + 0.49, nrows - ri - 0.51,
                    text,
                    ha="center", va="center",
                    fontsize=6.5, fontweight="bold", color="white",
                )

    ax.set_xlim(0, ncols)
    ax.set_ylim(0, nrows)
    ax.set_xticks([i + 0.49 for i in range(ncols)])
    ax.set_xticklabels(PLATE_COLS, fontsize=9)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_yticks([nrows - i - 0.51 for i in range(nrows)])
    ax.set_yticklabels(PLATE_ROWS, fontsize=9)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("Plate Layout", fontsize=13, fontweight="bold", pad=22)

    legend_handles = [mpatches.Patch(facecolor=palette[g], label=g)
                      for g in groups if g in palette]
    legend_handles.append(mpatches.Patch(facecolor="#888888", hatch="///", label="+ CytD treated"))
    ax.legend(handles=legend_handles, loc="upper left",
              bbox_to_anchor=(1.02, 1.0), fontsize=8,
              title="Groups", title_fontsize=9)

    fig.subplots_adjust(left=0.06, right=0.78, top=0.88, bottom=0.05)
    return fig


# ---------------------------------------------------------------------------
# Single-file category page (violin + strip)
# ---------------------------------------------------------------------------

def plot_category_page(
    df: pd.DataFrame,
    columns: list,
    title: str,
    groups: list,
    palette: dict,
    bacterial_page: bool = False,
) -> plt.Figure:
    """One subplot per metric column — violin + strip, grouped by treatment."""
    if bacterial_page:
        plot_df = df[df["has_bacteria"].notna()].copy()
        if plot_df.empty:
            plot_df = df.copy()
    else:
        plot_df = df.copy()

    available = [c for c in columns if c in plot_df.columns]
    if not available:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, f"No columns available for: {title}",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        return fig

    n = len(available)
    ncols = min(n, 4)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows), squeeze=False)
    fig.suptitle(title, fontsize=15, fontweight="bold", y=1.01)

    for idx, col in enumerate(available):
        ax = axes[idx // ncols][idx % ncols]
        col_data = plot_df[["group", col]].dropna(subset=[col])
        col_data = col_data[col_data["group"].isin(groups)]

        if col_data.empty:
            ax.set_visible(False)
            continue

        grp_order = [g for g in groups if g in col_data["group"].unique()]

        sns.violinplot(
            data=col_data, x="group", y=col,
            hue="group", order=grp_order, hue_order=grp_order,
            palette=palette, inner=None, linewidth=0.8, cut=0,
            legend=False, ax=ax,
        )
        sns.stripplot(
            data=col_data, x="group", y=col,
            order=grp_order, color="black", alpha=0.25, size=2.0,
            jitter=True, ax=ax,
        )

        ax.set_xlabel("")
        ax.set_ylabel(_nice_label(col), fontsize=9)
        ax.set_title(_nice_label(col), fontsize=10, fontweight="bold")
        ax.tick_params(axis="x", labelrotation=35, labelsize=8)
        ax.tick_params(axis="y", labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        for i, g in enumerate(grp_order):
            n_g = col_data[col_data["group"] == g].shape[0]
            ax.text(i, ax.get_ylim()[0], f"n={n_g}",
                    ha="center", va="bottom", fontsize=6.5, color="grey")

    for idx in range(len(available), nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Paired category page (grouped boxplots)
# ---------------------------------------------------------------------------

def plot_paired_category_page(
    df: pd.DataFrame,
    columns: list,
    title: str,
    comparison_label: str,
    groups: list,
    cond_order: list,
    cond_palette: dict,
    bacterial_page: bool = False,
) -> plt.Figure:
    """Violin + box plots: x = bacterial group, hue = condition. One subplot per metric."""
    if bacterial_page:
        plot_df = df[df["has_bacteria"].notna()].copy()
        if plot_df.empty:
            plot_df = df.copy()
    else:
        plot_df = df.copy()

    available = [c for c in columns if c in plot_df.columns]
    if not available:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, f"No columns available for: {title}",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title(f"{title}  [{comparison_label}]")
        return fig

    n = len(available)
    ncols = min(n, 4)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 5 * nrows), squeeze=False)
    fig.suptitle(f"{title}  [{comparison_label}]", fontsize=15, fontweight="bold", y=1.01)

    _BOX_WIDTH = 0.6
    n_conds = len(cond_order)

    for idx, col in enumerate(available):
        ax = axes[idx // ncols][idx % ncols]
        col_data = plot_df[["group", "condition", col]].dropna(subset=[col])
        col_data = col_data[col_data["group"].isin(groups)]

        if col_data.empty:
            ax.set_visible(False)
            continue

        grp_order = [g for g in groups if g in col_data["group"].unique()]

        # Violin (distribution shape, semi-transparent)
        sns.violinplot(
            data=col_data, x="group", y=col,
            hue="condition",
            order=grp_order,
            hue_order=cond_order,
            palette=cond_palette,
            width=_BOX_WIDTH,
            inner=None,
            cut=0,
            linewidth=0.5,
            alpha=0.38,
            dodge=True,
            legend=False,
            ax=ax,
        )

        # Box overlaid inside the violin (narrower, no outlier dots)
        sns.boxplot(
            data=col_data, x="group", y=col,
            hue="condition",
            order=grp_order,
            hue_order=cond_order,
            palette=cond_palette,
            width=_BOX_WIDTH * 0.28,
            gap=0.1,
            linewidth=0.9,
            fliersize=0,
            dodge=True,
            legend=False,
            ax=ax,
        )

        ax.set_xlabel("")
        ax.set_ylabel(_nice_label(col), fontsize=9)
        ax.set_title(_nice_label(col), fontsize=10, fontweight="bold")
        ax.tick_params(axis="x", labelrotation=35, labelsize=8)
        ax.tick_params(axis="y", labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # n= annotations aligned under each dodged box
        ymin = ax.get_ylim()[0]
        cond_color_list = list(cond_palette.values())
        for gi, g in enumerate(grp_order):
            for ci, cond in enumerate(cond_order):
                x_pos = gi + (ci - (n_conds - 1) / 2) * (_BOX_WIDTH / n_conds)
                n_val = col_data[
                    (col_data["group"] == g) & (col_data["condition"] == cond)
                ].shape[0]
                ax.text(x_pos, ymin, f"n={n_val}",
                        ha="center", va="bottom", fontsize=6,
                        color=cond_color_list[ci])

    for idx in range(len(available), nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    # Single figure-level condition legend
    legend_handles = [mpatches.Patch(facecolor=cond_palette[c], label=c) for c in cond_order]
    fig.legend(handles=legend_handles, loc="upper right",
               bbox_to_anchor=(1.0, 1.0), fontsize=9, title="Condition")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Single-file summary heatmap
# ---------------------------------------------------------------------------

def plot_summary_page(df: pd.DataFrame, groups: list, palette: dict) -> plt.Figure:
    """Z-scored group-mean heatmap across all numeric metrics."""
    _SKIP = {"experiment", "image_name", "cell_id", "group",
             "processing_mode", "measurement_source", "radial_method",
             "has_bacteria", "condition"}
    numeric_cols = [c for c in df.columns
                    if c not in _SKIP and pd.api.types.is_numeric_dtype(df[c])]

    group_means = df.groupby("group")[numeric_cols].mean()
    group_means = group_means.loc[[g for g in groups if g in group_means.index]]
    z = (group_means - group_means.mean()) / (group_means.std() + 1e-9)
    z = z.dropna(axis=1)

    fig, ax = plt.subplots(figsize=(max(14, len(z.columns) * 0.35), max(5, len(groups) * 0.8)))
    sns.heatmap(z, ax=ax, cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5,
                linewidths=0.3, linecolor="white",
                xticklabels=True, yticklabels=True,
                cbar_kws={"label": "Z-score (group mean)", "shrink": 0.5})
    ax.set_title("Group Mean Z-Score Heatmap — All Metrics", fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=60, labelsize=7)
    ax.tick_params(axis="y", rotation=0, labelsize=9)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Paired summary heatmap  (z-scores + delta)
# ---------------------------------------------------------------------------

def plot_paired_summary_page(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    label1: str,
    label2: str,
    groups: list,
    comparison_label: str,
    metric_cols: list,
) -> plt.Figure:
    """Two-panel heatmap: pooled z-scores (top) and delta z-score (bottom)."""
    shared_groups = [g for g in groups
                     if g in df1["group"].unique() and g in df2["group"].unique()]
    valid_cols = [c for c in metric_cols
                  if c in df1.columns and c in df2.columns
                  and pd.api.types.is_numeric_dtype(df1[c])]

    m1 = df1.groupby("group")[valid_cols].mean().reindex(shared_groups)
    m2 = df2.groupby("group")[valid_cols].mean().reindex(shared_groups)

    pooled = pd.concat([m1, m2], keys=[label1, label2])
    pooled_z = (pooled - pooled.mean()) / (pooled.std() + 1e-9)
    pooled_z = pooled_z.dropna(axis=1)
    z1 = pooled_z.loc[label1]
    z2 = pooled_z.loc[label2]
    delta = z2 - z1

    n_grp = len(shared_groups)
    fig_w = max(16, len(pooled_z.columns) * 0.45)
    fig, axes = plt.subplots(2, 1, figsize=(fig_w, 12),
                             gridspec_kw={"height_ratios": [2, 1], "hspace": 0.55})
    fig.suptitle(f"Summary Heatmaps — {comparison_label}", fontsize=14, fontweight="bold")

    sns.heatmap(pooled_z, ax=axes[0], cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5,
                linewidths=0.3, linecolor="white",
                xticklabels=True, yticklabels=True,
                cbar_kws={"label": "Z-score (group mean)", "shrink": 0.4})
    axes[0].set_title(
        f"Group Mean Z-Scores — {label1} (rows 1–{n_grp}) / {label2} (rows {n_grp + 1}–{2 * n_grp})",
        fontsize=9,
    )
    axes[0].tick_params(axis="x", rotation=60, labelsize=7)
    axes[0].tick_params(axis="y", rotation=0, labelsize=8)

    vabs = min(float(delta.abs().max().max()), 3.0)
    if np.isnan(vabs) or vabs == 0:
        vabs = 1.0
    sns.heatmap(delta, ax=axes[1], cmap="RdBu_r", center=0, vmin=-vabs, vmax=vabs,
                linewidths=0.3, linecolor="white",
                xticklabels=True, yticklabels=True,
                cbar_kws={"label": f"\u0394 Z-score ({label2} \u2212 {label1})", "shrink": 0.4})
    axes[1].set_title(
        f"Delta Z-Score ({label2} \u2212 {label1}) \u2014 positive = increased in {label2}",
        fontsize=9,
    )
    axes[1].tick_params(axis="x", rotation=60, labelsize=7)
    axes[1].tick_params(axis="y", rotation=0, labelsize=8)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Analytical pages
# ---------------------------------------------------------------------------

def plot_ranked_changes_page(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    label1: str,
    label2: str,
    metric_cols: list,
    groups: list,
    comparison_label: str,
) -> plt.Figure:
    """Horizontal bar chart of metrics ranked by |mean delta z-score| + per-group heatmap."""
    shared_groups = [g for g in groups
                     if g in df1["group"].unique() and g in df2["group"].unique()]
    valid_cols = [c for c in metric_cols
                  if c in df1.columns and c in df2.columns
                  and pd.api.types.is_numeric_dtype(df1[c])]

    m1 = df1.groupby("group")[valid_cols].mean().reindex(shared_groups)
    m2 = df2.groupby("group")[valid_cols].mean().reindex(shared_groups)
    pooled = pd.concat([m1, m2])
    pooled_z = (pooled - pooled.mean()) / (pooled.std() + 1e-9)
    z1 = pooled_z.iloc[:len(shared_groups)]
    z2 = pooled_z.iloc[len(shared_groups):]
    delta_df = (z2 - z1).copy()
    delta_df.index = shared_groups

    mean_delta = delta_df.mean(axis=0).dropna()
    mean_delta_sorted = mean_delta.reindex(mean_delta.abs().sort_values(ascending=False).index)

    n_metrics = len(mean_delta_sorted)
    fig, axes = plt.subplots(1, 2, figsize=(18, max(7, n_metrics * 0.32)))
    fig.suptitle(f"Ranked Metric Changes \u2014 {comparison_label}", fontsize=13, fontweight="bold")

    # Left: pooled bar chart
    ax = axes[0]
    bar_colors = ["#e74c3c" if v > 0 else "#3498db" for v in mean_delta_sorted.values]
    y_pos = list(range(n_metrics))
    ax.barh(y_pos, mean_delta_sorted.values, color=bar_colors, height=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([_nice_label(c) for c in mean_delta_sorted.index], fontsize=7)
    ax.invert_yaxis()
    ax.axvline(0, color="black", linewidth=0.8, zorder=3)
    ax.set_xlabel(f"Mean \u0394 Z-score ({label2} \u2212 {label1})", fontsize=9)
    ax.set_title("All Groups Pooled", fontsize=10, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Right: per-group delta heatmap (same metric order)
    delta_plot = delta_df[mean_delta_sorted.index]
    vabs = min(float(delta_plot.abs().max().max()), 3.0)
    if np.isnan(vabs) or vabs == 0:
        vabs = 1.0
    sns.heatmap(
        delta_plot.T,
        ax=axes[1],
        cmap="RdBu_r",
        center=0, vmin=-vabs, vmax=vabs,
        linewidths=0.3, linecolor="white",
        xticklabels=True, yticklabels=True,
        cbar_kws={"label": "\u0394 Z-score", "shrink": 0.5},
    )
    axes[1].set_title("Per-Group Breakdown", fontsize=10, fontweight="bold")
    axes[1].tick_params(axis="x", rotation=30, labelsize=8)
    axes[1].tick_params(axis="y", rotation=0, labelsize=7)

    fig.tight_layout()
    return fig


def plot_pca_page(
    df: pd.DataFrame,
    groups: list,
    palette: dict,
    comparison_label: str = "",
    cond_order: list = None,
    cond_palette: dict = None,
    metric_cols: list = None,
):
    """PCA scatter of metric columns. Marker shape = condition in paired mode.
    Returns None (and prints a warning) if scikit-learn is not installed."""
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        print("  [skip] PCA page: scikit-learn not installed (pip install scikit-learn)")
        return None

    if metric_cols is None:
        metric_cols = _all_metric_cols(df)

    cols_present = [c for c in metric_cols
                    if c in df.columns and pd.api.types.is_numeric_dtype(df[c])]
    extra_cols = ["group"] + (["condition"] if "condition" in df.columns else [])
    plot_df = df[extra_cols + cols_present].dropna().reset_index(drop=True)

    if len(plot_df) < 3 or len(cols_present) < 2:
        return None

    X_scaled = StandardScaler().fit_transform(plot_df[cols_present].values)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(X_scaled)
    var_exp = pca.explained_variance_ratio_ * 100

    _markers = ["o", "s", "D", "^", "v"]
    fig, ax = plt.subplots(figsize=(9, 7))

    if cond_order and "condition" in plot_df.columns:
        for ci, cond in enumerate(cond_order):
            for g in groups:
                mask = (plot_df["group"] == g) & (plot_df["condition"] == cond)
                if mask.sum() == 0:
                    continue
                ax.scatter(
                    coords[mask.values, 0], coords[mask.values, 1],
                    color=palette.get(g, "#888888"),
                    marker=_markers[ci % len(_markers)],
                    alpha=0.45, s=18, linewidths=0,
                )

        color_handles = [mpatches.Patch(facecolor=palette.get(g, "#888"), label=g)
                         for g in groups if g in plot_df["group"].values]
        from matplotlib.lines import Line2D
        shape_handles = [
            Line2D([0], [0], marker=_markers[ci % len(_markers)], color="grey",
                   linestyle="None", markersize=8, label=cond)
            for ci, cond in enumerate(cond_order)
        ]
        leg1 = ax.legend(handles=color_handles, loc="upper right", fontsize=8, title="Group")
        ax.add_artist(leg1)
        ax.legend(handles=shape_handles, loc="lower right", fontsize=8, title="Condition")
    else:
        for g in groups:
            mask = plot_df["group"] == g
            if mask.sum() == 0:
                continue
            ax.scatter(coords[mask.values, 0], coords[mask.values, 1],
                       color=palette.get(g, "#888888"),
                       alpha=0.45, s=18, linewidths=0, label=g)
        ax.legend(loc="upper right", fontsize=8, title="Group")

    ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)", fontsize=11)
    ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)", fontsize=11)
    ax.set_title(
        f"PCA Overview \u2014 {comparison_label}" if comparison_label else "PCA Overview",
        fontsize=12, fontweight="bold",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def plot_significance_page(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    label1: str,
    label2: str,
    groups: list,
    metric_cols: list,
    comparison_label: str,
) -> plt.Figure:
    """Mann-Whitney U heatmap: rows=metrics, cols=groups, colour = -log10(p).
    Raises ImportError if scipy is not installed."""
    from scipy.stats import mannwhitneyu

    shared_groups = [g for g in groups
                     if g in df1["group"].unique() and g in df2["group"].unique()]
    valid_cols = [c for c in metric_cols
                  if c in df1.columns and pd.api.types.is_numeric_dtype(df1[c])]

    pval_matrix = np.ones((len(valid_cols), len(shared_groups)))
    for mi, col in enumerate(valid_cols):
        for gi, g in enumerate(shared_groups):
            a = df1[df1["group"] == g][col].dropna().values
            b = df2[df2["group"] == g][col].dropna().values
            if len(a) >= 3 and len(b) >= 3:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    _, p = mannwhitneyu(a, b, alternative="two-sided")
                pval_matrix[mi, gi] = p

    log_p = -np.log10(pval_matrix + 1e-300)
    log_p_df = pd.DataFrame(
        log_p,
        index=[_nice_label(c) for c in valid_cols],
        columns=shared_groups,
    )

    fig_h = max(7, len(valid_cols) * 0.38)
    fig_w = max(8, len(shared_groups) * 1.3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        log_p_df, ax=ax,
        cmap="YlOrRd", vmin=0, vmax=5,
        linewidths=0.3, linecolor="white",
        xticklabels=True, yticklabels=True,
        cbar_kws={"label": "\u2212log\u2081\u2080(p-value)", "shrink": 0.5},
    )

    for mi in range(len(valid_cols)):
        for gi in range(len(shared_groups)):
            p = pval_matrix[mi, gi]
            star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            if star:
                ax.text(gi + 0.5, mi + 0.5, star,
                        ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")

    ax.set_title(f"Mann-Whitney U \u2014 {comparison_label}", fontsize=12, fontweight="bold")
    ax.tick_params(axis="x", rotation=30, labelsize=9)
    ax.tick_params(axis="y", rotation=0, labelsize=7)
    fig.tight_layout()
    return fig


def plot_correlation_page(
    df: pd.DataFrame,
    metric_cols: list,
    title: str = "",
    df2: pd.DataFrame = None,
    label1: str = "",
    label2: str = "",
) -> plt.Figure:
    """Pearson correlation heatmap(s). Side-by-side in paired mode."""
    valid_cols = [c for c in metric_cols
                  if c in df.columns and pd.api.types.is_numeric_dtype(df[c])]

    def _corr_heatmap(data, ax, subtitle):
        corr = data[valid_cols].corr(method="pearson")
        nice = [_nice_label(c) for c in corr.columns]
        corr.columns = nice
        corr.index = nice
        sns.heatmap(corr, ax=ax, cmap="coolwarm", center=0, vmin=-1, vmax=1,
                    linewidths=0.2, linecolor="white",
                    xticklabels=True, yticklabels=True,
                    cbar_kws={"label": "Pearson r", "shrink": 0.5})
        ax.set_title(subtitle, fontsize=10, fontweight="bold")
        ax.tick_params(axis="x", rotation=60, labelsize=6)
        ax.tick_params(axis="y", rotation=0, labelsize=6)

    n = len(valid_cols)
    if df2 is not None:
        fig_w = max(18, n * 0.55)
        fig_h = max(8, n * 0.45)
        fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h))
        _corr_heatmap(df, axes[0], f"Correlations \u2014 {label1}")
        _corr_heatmap(df2, axes[1], f"Correlations \u2014 {label2}")
        fig.suptitle(f"Metric Correlations \u2014 {title}", fontsize=13, fontweight="bold")
    else:
        fig_sz = max(10, n * 0.45)
        fig, ax = plt.subplots(figsize=(fig_sz, fig_sz))
        _corr_heatmap(df, ax, f"Metric Correlations \u2014 {title}")
        fig.suptitle(f"Metric Correlations \u2014 {title}", fontsize=13, fontweight="bold")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Combined experiment page  (two experiments side by side, shared y-axis per metric)
# ---------------------------------------------------------------------------

def plot_combined_category_page(
    exp_list: list,
    columns: list,
    title: str,
    groups: list,
    bacterial_page: bool = False,
) -> plt.Figure:
    """
    Side-by-side comparison of N experiments for one metric category.

    Layout: n_metrics rows x n_experiments columns.
    Y-axis is shared across columns per row so effect sizes are directly comparable.

    Each entry in exp_list must be a dict with keys:
        combined_df, cond_order, cond_palette, comparison_label
    """
    ref_df = exp_list[0]["combined_df"]
    available = [c for c in columns if c in ref_df.columns]
    if not available:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, f"No columns available for: {title}",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title(f"{title} — Combined")
        return fig

    n_metrics = len(available)
    n_exp = len(exp_list)
    _BOX_WIDTH = 0.6

    fig, axes = plt.subplots(
        n_metrics, n_exp,
        figsize=(7 * n_exp, 4.5 * n_metrics),
        squeeze=False,
        sharey="row",   # same y-scale per metric row — enables direct effect-size comparison
    )
    fig.suptitle(f"{title} — Combined View", fontsize=15, fontweight="bold", y=1.01)

    for mi, col in enumerate(available):
        for ei, exp in enumerate(exp_list):
            ax = axes[mi][ei]
            df = exp["combined_df"]
            cond_order = exp["cond_order"]
            cond_palette = exp["cond_palette"]

            if bacterial_page:
                df = df[df["has_bacteria"].notna()]
                if df.empty:
                    df = exp["combined_df"]

            col_data = df[["group", "condition", col]].dropna(subset=[col])
            col_data = col_data[col_data["group"].isin(groups)]

            if col_data.empty:
                ax.set_visible(False)
                continue

            grp_order = [g for g in groups if g in col_data["group"].unique()]
            n_conds = len(cond_order)

            # Violin
            sns.violinplot(
                data=col_data, x="group", y=col,
                hue="condition", order=grp_order, hue_order=cond_order,
                palette=cond_palette, width=_BOX_WIDTH,
                inner=None, cut=0, linewidth=0.5, alpha=0.38,
                dodge=True, legend=False, ax=ax,
            )
            # Box overlay
            sns.boxplot(
                data=col_data, x="group", y=col,
                hue="condition", order=grp_order, hue_order=cond_order,
                palette=cond_palette, width=_BOX_WIDTH * 0.28,
                gap=0.1, linewidth=0.9, fliersize=0,
                dodge=True, legend=False, ax=ax,
            )

            # Column header on top row only
            if mi == 0:
                ax.set_title(exp["comparison_label"], fontsize=10,
                             fontweight="bold", color="dimgrey", pad=8)

            # Y-label on leftmost column only
            if ei == 0:
                ax.set_ylabel(_nice_label(col), fontsize=9)
            else:
                ax.set_ylabel("")
                ax.tick_params(labelleft=False)

            ax.set_xlabel("")
            ax.tick_params(axis="x", labelrotation=35, labelsize=8)
            ax.tick_params(axis="y", labelsize=8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            # n= annotations
            ymin = ax.get_ylim()[0]
            cond_colors = list(cond_palette.values())
            for gi, g in enumerate(grp_order):
                for ci, cond in enumerate(cond_order):
                    x_pos = gi + (ci - (n_conds - 1) / 2) * (_BOX_WIDTH / n_conds)
                    n_val = col_data[
                        (col_data["group"] == g) & (col_data["condition"] == cond)
                    ].shape[0]
                    ax.text(x_pos, ymin, f"n={n_val}",
                            ha="center", va="bottom", fontsize=5.5,
                            color=cond_colors[ci])

    # Per-column condition legend on the bottom row
    for ei, exp in enumerate(exp_list):
        handles = [mpatches.Patch(facecolor=exp["cond_palette"][c], label=c)
                   for c in exp["cond_order"]]
        axes[-1][ei].legend(handles=handles, loc="upper right", fontsize=8, title="Condition")

    fig.tight_layout()
    return fig


def run_combined_experiment(cfg: dict):
    """Load all pairs in cfg['pairs'], render a combined side-by-side comparison PDF."""
    output_name = cfg["output_name"]
    pairs = cfg["pairs"]
    print(f"\nLoading (combined): {output_name}")

    exp_list = []
    for i, pair in enumerate(pairs):
        df1 = pd.read_csv(pair["csv1"])
        df2 = pd.read_csv(pair["csv2"])
        print(f"  Pair {i + 1} — {pair['label1']}: {len(df1)} rows, "
              f"{pair['label2']}: {len(df2)} rows")

        combined_preview = pd.concat([df1, df2], ignore_index=True)
        metric_cols = _all_metric_cols(combined_preview)
        metric_cols = _drop_constant_cols(combined_preview, metric_cols)

        df1 = _remove_outliers_iqr(df1, metric_cols)
        df2 = _remove_outliers_iqr(df2, metric_cols)
        df1["condition"] = pair["label1"]
        df2["condition"] = pair["label2"]
        combined = pd.concat([df1, df2], ignore_index=True)

        exp_list.append({
            "combined_df": combined,
            "cond_order": [pair["label1"], pair["label2"]],
            "cond_palette": {pair["label1"]: _COND_COLORS[i * 2 % len(_COND_COLORS)],
                             pair["label2"]: _COND_COLORS[(i * 2 + 1) % len(_COND_COLORS)]},
            "comparison_label": pair.get("comparison_label",
                                         f"{pair['label1']} vs {pair['label2']}"),
        })

    groups = _group_order(exp_list[0]["combined_df"])
    out_pdf = pathlib.Path(pairs[0]["csv1"]).parent / f"{output_name}_plots.pdf"

    with pdf_backend.PdfPages(out_pdf) as pdf:
        for category, columns in METRIC_GROUPS.items():
            # Only include columns present in every experiment
            cat_cols = [c for c in columns
                        if all(c in exp["combined_df"].columns for exp in exp_list)]
            if not cat_cols:
                continue
            print(f"  Rendering: {category}")
            is_bact = (category == "Bacterial Burden")
            fig = plot_combined_category_page(
                exp_list, cat_cols, category, groups, bacterial_page=is_bact
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        meta = pdf.infodict()
        meta["Title"] = f"TopoActin V3 \u2014 Combined: {output_name}"
        meta["Subject"] = "Combined cytoskeleton analysis metrics"

    print(f"  Saved: {out_pdf}")


# ---------------------------------------------------------------------------
# Experiment runners
# ---------------------------------------------------------------------------

def run_single_experiment(csv: pathlib.Path, output_name: str = None):
    """Load one CSV, apply QC, and render a single-file analysis PDF."""
    print(f"\nLoading (single): {csv.name}")
    df = pd.read_csv(csv)
    print(f"  {len(df)} rows, {df['group'].nunique()} groups: {list(df['group'].unique())}")

    all_metric_cols = _all_metric_cols(df)
    all_metric_cols = _drop_constant_cols(df, all_metric_cols)
    df = _remove_outliers_iqr(df, all_metric_cols)

    groups = _group_order(df)
    palette = _make_palette(groups)
    stem = output_name if output_name else csv.stem
    out_pdf = csv.parent / f"{stem}_plots.pdf"

    with pdf_backend.PdfPages(out_pdf) as pdf:
        if SHOW_PLATE_MAP:
            print("  Rendering: Plate Map")
            fig = plot_plate_map_page(groups, palette)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        print("  Rendering: Summary Heatmap")
        fig = plot_summary_page(df, groups, palette)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        for category, columns in METRIC_GROUPS.items():
            cat_cols = [c for c in columns if c in all_metric_cols]
            if not cat_cols:
                continue
            print(f"  Rendering: {category}")
            is_bact = (category == "Bacterial Burden")
            fig = plot_category_page(df, cat_cols, category, groups, palette,
                                     bacterial_page=is_bact)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        print("  Rendering: PCA Overview")
        fig = plot_pca_page(df, groups, palette, comparison_label=stem,
                            metric_cols=all_metric_cols)
        if fig is not None:
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        print("  Rendering: Metric Correlations")
        fig = plot_correlation_page(df, all_metric_cols, title=stem)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        meta = pdf.infodict()
        meta["Title"] = f"TopoActin V3 \u2014 {stem}"
        meta["Subject"] = "Cytoskeleton analysis metrics"

    print(f"  Saved: {out_pdf}")


def run_paired_experiment(
    csv1: pathlib.Path,
    csv2: pathlib.Path,
    label1: str,
    label2: str,
    output_name: str,
):
    """Load two CSVs, combine with condition labels, and render a paired analysis PDF."""
    comparison_label = f"{label1} vs {label2}"
    print(f"\nLoading (paired): {comparison_label}")
    df1 = pd.read_csv(csv1)
    df2 = pd.read_csv(csv2)
    print(f"  {label1}: {len(df1)} rows \u2014 {csv1.name}")
    print(f"  {label2}: {len(df2)} rows \u2014 {csv2.name}")

    # Determine effective metric columns across both files
    combined_preview = pd.concat([df1, df2], ignore_index=True)
    all_metric_cols = _all_metric_cols(combined_preview)
    all_metric_cols = _drop_constant_cols(combined_preview, all_metric_cols)

    # Outlier removal applied to each file independently
    df1 = _remove_outliers_iqr(df1, all_metric_cols)
    df2 = _remove_outliers_iqr(df2, all_metric_cols)

    # Add condition labels after cleaning
    df1["condition"] = label1
    df2["condition"] = label2
    combined = pd.concat([df1, df2], ignore_index=True)

    groups = _group_order(combined)
    palette = _make_palette(groups)
    cond_order = [label1, label2]
    cond_palette = {label1: _COND_COLORS[0], label2: _COND_COLORS[1]}
    out_pdf = csv1.parent / f"{output_name}_plots.pdf"

    with pdf_backend.PdfPages(out_pdf) as pdf:
        if SHOW_PLATE_MAP:
            print("  Rendering: Plate Map")
            fig = plot_plate_map_page(groups, palette)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        print("  Rendering: Paired Summary Heatmap")
        fig = plot_paired_summary_page(
            df1, df2, label1, label2, groups, comparison_label, all_metric_cols
        )
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        for category, columns in METRIC_GROUPS.items():
            cat_cols = [c for c in columns if c in all_metric_cols]
            if not cat_cols:
                continue
            print(f"  Rendering: {category}")
            is_bact = (category == "Bacterial Burden")
            fig = plot_paired_category_page(
                combined, cat_cols, category, comparison_label,
                groups, cond_order, cond_palette, bacterial_page=is_bact,
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        print("  Rendering: Ranked Metric Changes")
        fig = plot_ranked_changes_page(
            df1, df2, label1, label2, all_metric_cols, groups, comparison_label
        )
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        print("  Rendering: PCA Overview")
        fig = plot_pca_page(combined, groups, palette, comparison_label,
                            cond_order, cond_palette, all_metric_cols)
        if fig is not None:
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        print("  Rendering: Statistical Significance (Mann-Whitney U)")
        try:
            fig = plot_significance_page(
                df1, df2, label1, label2, groups, all_metric_cols, comparison_label
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
        except ImportError:
            print("  [skip] Significance page: scipy not installed (pip install scipy)")

        print("  Rendering: Metric Correlations")
        fig = plot_correlation_page(
            df1, all_metric_cols, title=comparison_label,
            df2=df2, label1=label1, label2=label2,
        )
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        meta = pdf.infodict()
        meta["Title"] = f"TopoActin V3 \u2014 {comparison_label}"
        meta["Subject"] = "Paired cytoskeleton analysis metrics"

    print(f"  Saved: {out_pdf}")


# ---------------------------------------------------------------------------
# Main dispatcher
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) >= 2:
        # Backwards-compatible: single CSV path as argument
        csv_path = pathlib.Path(sys.argv[1])
        output_name = sys.argv[2] if len(sys.argv) >= 3 else None
        run_single_experiment(csv_path, output_name)
        return

    if not EXPERIMENT_CONFIGS:
        print("No EXPERIMENT_CONFIGS defined and no CSV path provided.")
        print("Usage: python plot_results.py <path_to_csv>")
        sys.exit(1)

    for cfg in EXPERIMENT_CONFIGS:
        mode = cfg.get("mode", "single")
        if mode == "single":
            run_single_experiment(
                pathlib.Path(cfg["csv"]),
                cfg.get("output_name"),
            )
        elif mode == "paired":
            run_paired_experiment(
                pathlib.Path(cfg["csv1"]),
                pathlib.Path(cfg["csv2"]),
                cfg["label1"],
                cfg["label2"],
                cfg["output_name"],
            )
        elif mode == "combined":
            run_combined_experiment(cfg)
        else:
            print(f"Unknown mode '{mode}' in EXPERIMENT_CONFIGS — skipping.")

    print("\nAll experiments complete.")


if __name__ == "__main__":
    main()
