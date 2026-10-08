"""
analyze_columns.py — Identify redundant columns in TopoActin V3 pipeline CSV output.

Usage:
    python analyze_columns.py <path_to_csv>

Outputs:
    {csv_stem}_deduped.csv          — confirmed alias columns removed (55 → 49 cols)
    {csv_stem}_correlation_heatmap.png — Pearson correlation matrix of all numeric columns
    Terminal report with column groups and collinearity warnings
"""

import sys
import pathlib
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Confirmed alias columns — same computation, legacy duplicate name
# These are SAFE to auto-remove; the canonical name is kept.
# ---------------------------------------------------------------------------
ALIAS_REMOVALS = {
    "zone_1_density":        "inner_density",           # same perinuclear zone
    "zone_2_density":        "mid_density",              # same mid zone
    "zone_d_density":        "mid_density",              # second alias of mid zone
    "zone_3_density":        "outer_density",            # same cortical zone
    "zone_r_density":        "outer_density",            # second alias of cortical zone
    "intensity_capture_ratio": "fidelity",               # identical metric
}

# ---------------------------------------------------------------------------
# Columns to skip when computing correlations (non-numeric / identifiers)
# ---------------------------------------------------------------------------
SKIP_COLS = {
    "experiment", "image_name", "cell_id", "group",
    "processing_mode", "measurement_source", "radial_method", "has_bacteria",
}


def compute_correlation_matrix(df: pd.DataFrame) -> pd.DataFrame:
    num_cols = [c for c in df.columns if c not in SKIP_COLS
                and pd.api.types.is_numeric_dtype(df[c])]
    return df[num_cols].corr(method="pearson")


def find_high_corr_pairs(corr: pd.DataFrame, threshold: float = 0.90):
    pairs = []
    cols = corr.columns.tolist()
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            r = corr.iloc[i, j]
            if pd.notna(r) and abs(r) >= threshold:
                pairs.append((abs(r), cols[i], cols[j]))
    pairs.sort(reverse=True)
    return pairs


def save_heatmap(corr: pd.DataFrame, out_path: pathlib.Path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns

        size = max(12, len(corr) * 0.38)
        fig, ax = plt.subplots(figsize=(size, size * 0.85))
        mask = np.triu(np.ones_like(corr, dtype=bool))
        sns.heatmap(
            corr,
            mask=mask,
            ax=ax,
            cmap="RdBu_r",
            center=0,
            vmin=-1, vmax=1,
            linewidths=0.3,
            linecolor="white",
            xticklabels=True,
            yticklabels=True,
            cbar_kws={"label": "Pearson r", "shrink": 0.6},
            annot=False,
        )
        ax.set_title("Column Correlation Matrix — TopoActin V3 Metrics", fontsize=13, fontweight="bold")
        ax.tick_params(axis="x", rotation=60, labelsize=7)
        ax.tick_params(axis="y", rotation=0, labelsize=7)
        fig.tight_layout()
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return True
    except ImportError:
        return False


def main(csv_path: pathlib.Path):
    print(f"\nLoading: {csv_path}")
    df = pd.read_csv(csv_path)
    print(f"  {len(df)} rows × {len(df.columns)} columns\n")

    # ── 1. Alias removal ────────────────────────────────────────────────────
    present_aliases = {col: canon for col, canon in ALIAS_REMOVALS.items()
                       if col in df.columns}
    missing_canon   = [canon for canon in present_aliases.values()
                       if canon not in df.columns]

    print("=" * 70)
    print("GROUP A — Confirmed aliases (safe to remove)")
    print("=" * 70)
    if present_aliases:
        for col, canon in present_aliases.items():
            if canon in df.columns:
                # Verify they really are identical
                both = df[[col, canon]].dropna()
                if len(both) > 0:
                    r = both[col].corr(both[canon])
                    identical = np.allclose(both[col].values, both[canon].values, equal_nan=True)
                    status = "IDENTICAL" if identical else f"r={r:.4f}"
                else:
                    status = "no overlap to verify"
                print(f"  REMOVE  {col:<35}  (≡ {canon})   [{status}]")
            else:
                print(f"  SKIP    {col:<35}  (canonical '{canon}' not in CSV)")
    else:
        print("  None of the known alias columns found in this CSV.")

    df_deduped = df.drop(columns=[c for c in present_aliases if c in df.columns
                                  and present_aliases[c] in df.columns],
                         errors="ignore")
    removed = len(df.columns) - len(df_deduped.columns)
    print(f"\n  Removed {removed} alias columns: {len(df.columns)} → {len(df_deduped.columns)} columns")

    # ── 2. Correlation analysis ─────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("GROUP B — High-correlation pairs (|r| > 0.90) — flagged, NOT removed")
    print("=" * 70)
    corr = compute_correlation_matrix(df_deduped)
    high_pairs = find_high_corr_pairs(corr, threshold=0.90)

    # Filter out any remaining aliases that weren't removed (shouldn't happen)
    alias_set = set(ALIAS_REMOVALS.keys())
    high_pairs = [(r, c1, c2) for r, c1, c2 in high_pairs
                  if c1 not in alias_set and c2 not in alias_set]

    if high_pairs:
        print(f"  {'r':>7}  {'Column A':<38}  {'Column B'}")
        print(f"  {'-'*7}  {'-'*38}  {'-'*38}")
        for r, c1, c2 in high_pairs[:25]:
            flag = " ← PERFECT" if r > 0.9999 else ""
            print(f"  {r:7.4f}  {c1:<38}  {c2}{flag}")
        if len(high_pairs) > 25:
            print(f"  ... and {len(high_pairs) - 25} more pairs")
        print(f"\n  These pairs are flagged for your awareness.")
        print(f"  Do NOT remove without domain knowledge — they may diverge in other datasets.")
    else:
        print("  No high-correlation pairs found after alias removal.")

    # ── 3. Save deduplicated CSV ────────────────────────────────────────────
    out_csv = csv_path.with_name(csv_path.stem + "_deduped.csv")
    df_deduped.to_csv(out_csv, index=False)
    print(f"\nSaved deduplicated CSV: {out_csv}")

    # ── 4. Correlation heatmap ──────────────────────────────────────────────
    out_png = csv_path.with_name(csv_path.stem + "_correlation_heatmap.png")
    ok = save_heatmap(corr, out_png)
    if ok:
        print(f"Saved correlation heatmap:  {out_png}")
    else:
        print("Heatmap skipped (matplotlib/seaborn not available).")

    # ── 5. Summary ──────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Original columns:          {len(df.columns)}")
    print(f"  Alias columns removed:     {removed}")
    print(f"  Remaining columns:         {len(df_deduped.columns)}")
    print(f"  High-corr pairs flagged:   {len(high_pairs)}")
    print(f"\n  Retained columns:")
    for i, c in enumerate(df_deduped.columns):
        print(f"    {i:2d}: {c}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyze_columns.py <path_to_csv>")
        sys.exit(1)
    main(pathlib.Path(sys.argv[1]))
