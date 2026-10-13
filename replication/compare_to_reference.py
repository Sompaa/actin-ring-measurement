"""Compare a re-run with the original paper CSV, cell by cell and well by well.

    python replication/compare_to_reference.py NEW.csv REFERENCE.csv -o compare_out/

Only the fields present in NEW.csv are compared. cell_id is not stable between
runs, so cells are paired within each field by area and perimeter (Hungarian
matching; a pair is accepted when the two log-ratios add up to < 0.10).

Writes:
  cells_per_field.csv   cell count per field, new vs reference
  metric_agreement.csv  per metric: matched cells, median ratio new/ref,
                        correlation, median relative difference
  well_medians.csv      per well and metric: median new vs median reference
                        (all cells, not only matched ones), useful on whole plates
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

SKIP = {"experiment", "image_name", "cell_id", "group", "processing_mode", "radial_method", "measurement_source"}


def match_cells(ref, new, images, max_cost=0.10):
    pairs = []
    for im in images:
        r = ref[ref.image_name == im].reset_index(drop=True)
        n = new[new.image_name == im].reset_index(drop=True)
        if r.empty or n.empty:
            continue
        cost = (np.abs(np.log(r.cell_area.values[:, None] / n.cell_area.values[None]))
                + np.abs(np.log(r.cell_perimeter.values[:, None] / n.cell_perimeter.values[None])))
        ri, ni = linear_sum_assignment(cost)
        pairs += [(r.iloc[a], n.iloc[b]) for a, b in zip(ri, ni) if cost[a, b] < max_cost]
    return pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("new_csv")
    ap.add_argument("reference_csv")
    ap.add_argument("-o", "--out", default="compare_out")
    args = ap.parse_args()

    new = pd.read_csv(args.new_csv, low_memory=False)
    ref = pd.read_csv(args.reference_csv, low_memory=False)
    images = sorted(set(new.image_name) & set(ref.image_name))
    only_new = sorted(set(new.image_name) - set(ref.image_name))
    ref = ref[ref.image_name.isin(images)]
    new_c = new[new.image_name.isin(images)]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    counts = pd.DataFrame({"reference": ref.groupby("image_name").size(),
                           "new": new_c.groupby("image_name").size()}).fillna(0).astype(int)
    counts.to_csv(out / "cells_per_field.csv")

    metrics = [c for c in ref.columns if c not in SKIP and c in new.columns
               and pd.api.types.is_numeric_dtype(ref[c])]
    pairs = match_cells(ref, new_c, images)
    rows = []
    for c in metrics:
        r = np.array([p[0][c] for p in pairs], float)
        n = np.array([p[1][c] for p in pairs], float)
        ok = np.isfinite(r) & np.isfinite(n) & (r != 0)
        if ok.sum() < 3:
            continue
        corr = np.corrcoef(r[ok], n[ok])[0, 1] if r[ok].std() > 0 and n[ok].std() > 0 else np.nan
        rows.append({"metric": c, "matched_cells": int(ok.sum()),
                     "median_ratio_new_over_ref": np.median(n[ok] / r[ok]),
                     "correlation": corr,
                     "median_rel_diff": np.median(np.abs(n[ok] - r[ok]) / np.abs(r[ok]))})
    agree = pd.DataFrame(rows, columns=["metric", "matched_cells", "median_ratio_new_over_ref",
                                        "correlation", "median_rel_diff"])
    agree.to_csv(out / "metric_agreement.csv", index=False)

    for d in (ref, new_c):
        d.loc[:, "well"] = d.image_name.str.split("_").str[0]
    wm = pd.concat({"reference": ref.groupby("well")[metrics].median(),
                    "new": new_c.groupby("well")[metrics].median()}, axis=1)
    wm.to_csv(out / "well_medians.csv")

    print(f"Fields compared: {len(images)}"
          + (f"  ({len(only_new)} field(s) in the new run have no reference rows)" if only_new else ""))
    print(f"Cells: reference {len(ref)}, new {len(new_c)}, matched by size/shape {len(pairs)}")
    focus = agree[agree.metric.isin(["cell_area", "cortical_intensity", "inner_intensity",
                                     "cortical_weighted_density", "skeleton_density"])]
    print(focus.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\nTables written to {out}/")


if __name__ == "__main__":
    main()
