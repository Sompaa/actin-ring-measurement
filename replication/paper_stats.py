"""Statistical analysis of one or more per-cell result tables.

    python replication/paper_stats.py RESULT.csv [MORE.csv ...] -o stats_out/

Each CSV is treated as one experiment and analysed on its own:

1. within each treatment x well, clip each value to that well's 5th-95th
   percentile (winsorizing: extreme cells are pulled in, not removed)
2. one median per feature per well (the well is the replicate, not the cell)
3. KpSen (ancestor) wells vs BasR wells: unpaired Mann-Whitney U test
4. effect = median of BasR well medians minus median of KpSen well medians,
   with a 95% interval from 2000 bootstrap resamples of the wells
5. Benjamini-Hochberg correction across all features of that experiment

Writes stats_<name>.csv per experiment and stats_headline.csv with
cortical_weighted_density and cortical_intensity side by side.

With two or more CSVs it also pools the experiments: one value per well, a linear mixed-effects
model  well_median ~ group + (1 | experiment)  for KpSen vs BasR, so that
plate-to-plate differences are absorbed by the experiment term, then
Benjamini-Hochberg across features -> stats_pooled_lmm.csv, and the per-well
medians behind the violin plots -> well_medians_pooled.csv.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HEADLINE = ["cortical_weighted_density", "cortical_intensity"]
META = {"experiment", "image_name", "cell_id", "group", "processing_mode", "radial_method",
        "measurement_source", "well", "Treatment"}


def bh(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(pv)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def well_medians(df, features, lo, hi):
    df = df.copy()
    df["well"] = df["image_name"].str.split("_").str[0]

    def clip(g):
        g = g.copy()
        for f in features:
            ql, qh = g[f].quantile([lo, hi])
            g[f] = g[f].clip(ql, qh)
        return g

    df = df.groupby(["group", "well"], group_keys=False)[features + ["group", "well"]].apply(clip)
    return df.groupby(["group", "well"])[features].median().reset_index()


def analyse(df, ref="KpSen", test="BasR", lo=0.05, hi=0.95, n_boot=2000, seed=0):
    features = [c for c in df.columns if c not in META and pd.api.types.is_numeric_dtype(df[c])]
    n_cells = len(df)
    wm = well_medians(df, features, lo, hi)
    a = wm[wm["group"] == ref]
    b = wm[wm["group"] == test]
    rng = np.random.default_rng(seed)
    rows = []
    for f in features:
        x, y = a[f].dropna().to_numpy(), b[f].dropna().to_numpy()
        if len(x) < 2 or len(y) < 2:
            continue
        p = mannwhitneyu(y, x, alternative="two-sided").pvalue
        diff = np.median(y) - np.median(x)
        boots = np.median(rng.choice(y, (n_boot, len(y))), axis=1) - np.median(rng.choice(x, (n_boot, len(x))), axis=1)
        lo_ci, hi_ci = np.percentile(boots, [2.5, 97.5])
        rows.append({"feature": f, f"{ref}_median": np.median(x), f"{test}_median": np.median(y),
                     "median_diff": diff, "ci95_low": lo_ci, "ci95_high": hi_ci,
                     f"n_wells_{ref}": len(x), f"n_wells_{test}": len(y), "p": p})
    res = pd.DataFrame(rows)
    res["q_BH"] = bh(res["p"])
    return res.sort_values("p").reset_index(drop=True), n_cells


def pooled_lmm(frames, ref="KpSen", test="BasR"):
    import warnings
    import statsmodels.formula.api as smf

    wms = []
    for name, df in frames:
        features = [c for c in df.columns if c not in META and pd.api.types.is_numeric_dtype(df[c])]
        wms.append(well_medians(df, features, 0.05, 0.95).assign(experiment=name))
    wm = pd.concat(wms, ignore_index=True)
    feats = [c for c in wm.columns if c not in {"group", "well", "experiment"}]
    sub = wm[wm.group.isin([ref, test])].copy()
    sub["is_test"] = (sub.group == test).astype(float)
    rows = []
    for f in feats:
        d = sub[["experiment", "is_test", f]].dropna().rename(columns={f: "y"})
        if d.y.std() == 0 or d.experiment.nunique() < 2:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = smf.mixedlm("y ~ is_test", d, groups=d["experiment"]).fit(reml=True)
        lo, hi = fit.conf_int().loc["is_test"]
        rows.append({"feature": f, f"{ref}_median": d.y[d.is_test == 0].median(),
                     f"{test}_median": d.y[d.is_test == 1].median(),
                     "lmm_effect": fit.params["is_test"], "ci95_low": lo, "ci95_high": hi,
                     "p": fit.pvalues["is_test"], "n_wells": len(d)})
    res = pd.DataFrame(rows)
    res["q_BH"] = bh(res["p"])
    return res.sort_values("p").reset_index(drop=True), wm


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="+")
    ap.add_argument("-o", "--out", default="stats_out")
    ap.add_argument("--names", nargs="*", help="short names, one per CSV (default: file stem)")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    names = args.names or [Path(c).stem[:40] for c in args.csv]
    head, frames = [], []
    for name, path in zip(names, args.csv):
        df = pd.read_csv(path, low_memory=False)
        frames.append((name, df))
        res, n = analyse(df)
        res.to_csv(out / f"stats_{name}.csv", index=False)
        print(f"\n== {name}: {n} cells")
        h = res[res["feature"].isin(HEADLINE)]
        print(h[["feature", "KpSen_median", "BasR_median", "median_diff", "ci95_low", "ci95_high", "p", "q_BH"]]
              .to_string(index=False, float_format=lambda v: f"{v:.4g}"))
        head.append(h.assign(experiment=name))
    pd.concat(head).to_csv(out / "stats_headline.csv", index=False)

    if len(frames) > 1:
        res, wm = pooled_lmm(frames)
        res.to_csv(out / "stats_pooled_lmm.csv", index=False)
        wm.to_csv(out / "well_medians_pooled.csv", index=False)
        print(f"\n== Pooled, mixed model with experiment as random effect ({len(frames)} experiments)")
        print(res[res.feature.isin(HEADLINE)][["feature", "KpSen_median", "BasR_median", "lmm_effect",
                                               "ci95_low", "ci95_high", "p", "q_BH"]]
              .to_string(index=False, float_format=lambda v: f"{v:.4g}"))
        print("\nWell-median ranges (violin plots, plot_figure_FG.py):")
        print(wm.groupby("group")[HEADLINE].agg(["min", "median", "max"]).round(1).to_string())
    print(f"\nTables written to {out}/")


if __name__ == "__main__":
    main()
