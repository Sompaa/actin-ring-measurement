# TopoActin V3 — Development Workbook
## February 2026 | A549 mScarlet Actin Cytoskeleton Pipeline

---

## Overview

This month was mostly about turning the pipeline from a pure actin topology tool into something that can actually speak to the biology — meaning we can now ask *where inside the cell are the bacteria, and what is the actin doing there?* On top of that we cleaned up a bunch of legacy metric duplication that had been creeping in, added a more principled image normalization step, and started building a validation framework using synthetic images. A lot of the work was incremental and somewhat reactive (bugs found during test runs, notebook UX issues), but the core additions are solid and ready for the first real batch analysis on the server.

---

## 1. Bacterial Burden Integration into the Actin Pipeline

### What we built

The main goal was to connect two previously separate analysis worlds: the actin topology pipeline (TopoActin V3) and the bacterial quantification module (`bacterial_burden/`). The result is a unified per-cell measurement that gives you, in the same CSV row, both the cytoskeletal fingerprint of a cell *and* how many bacteria are inside it and where they sit spatially.

The key outputs added to the pipeline CSV:

| Column | What it measures |
|---|---|
| `intracellular_bact_area` | Total bacteria pixel area inside the cell |
| `area_ratio` | Bacteria area / host cell area |
| `has_bacteria` | Boolean — any bacteria detected in this cell |
| `bact_density_perinuclear` | Bacteria density in the inner zone (0–33% of cell radius from nucleus) |
| `bact_density_middle` | Bacteria density in the mid zone (33–66%) |
| `bact_density_cortical` | Bacteria density in the outer/cortical zone (66–100%) |

The zone boundaries use exactly the same nucleus-anchored distance map as the actin topology analysis, so the bacterial and actin zones are directly comparable. This was intentional — you can now ask "does cortical bacterial density correlate with cortical actin density?" with spatially aligned measurements.

**Detection approach:** Sato ridge filter on the bacteria channel (optimized for rod-shaped bacteria: σ = 1.0–3.0 µm), followed by Otsu thresholding, per-cell percentile normalization (2nd–98th), and a minimum signal range QC check to avoid false positives in cells with essentially no bacteria channel signal.

**Channel mapping (corrected midway):** actin=0, nuclei=1, brightfield=2, bacteria=3 for 4-channel images. There was a channel swap bug found during the first test run — brightfield and bacteria were inverted — which was fixed across the config, image loader, and notebook.

**Batch mode:** Each experiment in a batch run gets its own output CSV with the full bacterial burden columns. Individual experiments can override the `bacterial_burden` block (e.g., disable bacteria for a 2-channel experiment in the same batch).

### Struggles

- **ImportError in Jupyter kernel**: After adding `compute_zone_bacterial_density` to `detection.py`, the notebook kept throwing `ImportError: cannot import name 'compute_zone_bacterial_density'`. Root cause: Jupyter had cached the old version of the module in the kernel. Fixed by adding `importlib.reload()` at the top of the bacterial burden preview cell. This won't be an issue in the main pipeline (fresh Python process each run).

- **NaN vs zero in zone columns**: Cells with bacteria signal below the minimum signal range threshold were getting `NaN` in all six zone columns instead of `0.0`. The difference matters for statistical analyses — `NaN` means "not measured", `0.0` means "measured and found nothing". Fixed by adding an empty (all-False) bacteria mask to the early-exit result dict, so `compute_zone_bacterial_density` still runs and returns zeros.

- **Channel swap discovery**: The first test run produced nonsensical bacterial burden numbers. Turned out brightfield and bacteria channels were transposed in the config template. Fixed in three places: config JSON, notebook image loader cell, and notebook config export cell.

---

## 2. Metric Redundancy Reduction and Image Normalization

### Column deduplication

The pipeline CSV had grown to 55 columns, with several being exact duplicates under different legacy names. These came from the pipeline's dual-path architecture having two separate naming conventions that both ended up in the output. Running a pairwise correlation analysis on the real dataset revealed:

**Exact aliases removed (55 → 49 columns):**
- `zone_1_density` → kept as `inner_density`
- `zone_2_density` + `zone_d_density` → kept as `mid_density`
- `zone_3_density` + `zone_r_density` → kept as `outer_density`
- `intensity_capture_ratio` → kept as `fidelity`

These had `r = 1.000` in all tested rows, confirming they're computationally identical.

**High-correlation pairs flagged but NOT removed:** `num_branches`/`num_junctions`/`total_network_intensity` all correlated near-perfectly in the A549 dataset (r ≈ 0.99), and `total_length` vs `weighted_total_length` were also near-identical. These are kept because they represent genuinely different measurements that could diverge in cell types with more variable actin intensity distributions. They're documented and flagged for collinearity awareness in downstream modeling.

The `analyze_columns.py` script can be run on any future CSV to auto-detect and remove confirmed aliases, and produces a correlation heatmap for manual inspection.

### Image-level CLAHE normalization

This was a conceptual clarification as much as an implementation task. The pipeline already had a CLAHE option in the skeletonization config (`use_clahe`), but it was applied *to the normalized Sato map* after the raw Sato values were already saved. This meant it had zero effect on the weighted density metrics (Path C) — only the binary skeleton topology changed.

The new image-level CLAHE is applied to the **full actin MIP** right after brightness/contrast adjustment, before any per-cell processing. This means:
- All cells in the same image share the same local contrast enhancement
- Both the binary skeleton (Path B) and intensity-weighted metrics (Path C) are affected
- A dim cell doesn't get artificially equalized to look like a bright cell — the enhancement is local, so spatial differences within the image are preserved

Config:
```json
"image_normalization": {
  "use_clahe": false,
  "clahe_clip_limit": 0.03,
  "clahe_tile_grid_size": [8, 8]
}
```

It's off by default. For the A549 mScarlet experiments, enabling it with clip_limit 0.02–0.04 is recommended given the typical vignetting and cell-to-cell brightness variation in these plates.

The interactive notebook was updated to reflect this: the CLAHE toggle in the brightness/contrast step now applies CLAHE *after* B/C (not instead of it), and a tile size slider was added. The `apply_clahe()` helper was also fixed — it was accepting a `tile_size` parameter but not actually passing it to scikit-image's `equalize_adapthist`.

### Why this matters

Without image-level normalization, the weighted density metrics (inner/mid/cortical actin density) are sensitive to absolute brightness differences between fields. A well with slightly different staining efficiency would produce systematically different density values even if the actin structure is the same. The combination of B/C preprocessing + optional CLAHE gives better within-experiment comparability, which is important for the multi-strain infection comparisons.

---

## 3. Synthetic Image Generation for Pipeline Validation

### Why synthetic data

The pipeline produces a rich set of metrics, but without ground truth it's hard to say whether the numbers are correct or just consistent. Synthetic images with *known* structure let you run the full pipeline and check whether the output makes biological and mathematical sense.

Five scenarios were created, each as a `[Z=5, C=4, Y=512, X=512]` uint16 TIFF matching the real pipeline input format:

| Scenario | Actin pattern | Bacteria | What it validates |
|---|---|---|---|
| `stress_fibers` | 10 parallel horizontal fibers | None | Tortuosity should be ~1.0; fibers are straight |
| `cortical_ring` | Thick peripheral ring, sparse center | None | Cortical density >> inner density |
| `dense_mesh` | 60 random-angle filaments | None | High branch count, moderate tortuosity |
| `infected_cell` | Cortical ring + radial fibers | 4 rods near nucleus | `bact_density_perinuclear` > 0, cortical ≈ 0 |
| `empty_control` | Near-background noise | None | Pipeline shouldn't crash on low-signal cells |

Each image includes realistic Gaussian noise (σ ≈ 200 counts on a ~500 count background), Z-defocus variation so the MIP projection works correctly, and a proper nucleus channel for zone mapping.

A `ground_truth.csv` is included with the expected metric directions for each scenario. The intention is to run the pipeline on this directory, compare the output CSV to ground truth, and quantify per-metric accuracy.

### Limitations to keep in mind

Synthetic images can validate the *direction* of metrics (cortical > inner for the ring scenario) but not the *absolute values* reliably — the pipeline applies several preprocessing steps (B/C adjustment, Sato filter, hysteresis thresholding) that make pixel-exact ground truth comparisons difficult. The more useful output from validation is whether the pipeline's *ranking* of scenarios matches expectations (e.g., dense_mesh should have more branches than stress_fibers) and whether pathological cases (empty_control) are handled gracefully.

---

## Struggles and Lessons Learned

Beyond the specific bugs listed above, a few broader patterns:

- **Jupyter notebook editing is painful at scale.** The `.ipynb` format means you can't use standard file editing tools — everything has to go through JSON manipulation or the `NotebookEdit` API. Several cells took multiple attempts because the exact source string didn't match due to hidden whitespace or previously unread modifications.

- **The dual-path architecture created metric confusion.** Path B (structure/binary) and Path C (intensity weighting) are computed separately and the columns ended up with overlapping names that didn't clearly indicate which path they came from. The alias cleanup helped, but a proper naming audit is still needed.

- **Windows + conda + Unicode = headaches.** The Windows terminal's CP1250 encoding broke the output of scripts that used Unicode characters (arrows, checkmarks, Greek letters). Scripts now use ASCII-only output for compatibility.

- **"The CLAHE is already implemented"** — except it wasn't doing what we thought. The config had `use_clahe` in the skeletonization section but it only affected binarization, not the metrics. This took careful tracing of the data flow to identify.

---

## Future Plans

### 4. Colocalization Study: Human Actin vs Bacterial Invasion Zones

The zone-based bacterial density metrics open up a colocalization analysis that wasn't possible before. The next step is to correlate, on a per-cell basis:
- Bacterial presence in each zone vs actin density in that zone
- Whether bacteria preferentially accumulate in high- or low-actin regions
- Whether infection history (time point, strain) shifts the bacterial distribution (perinuclear vs cortical)

Concretely, this means computing a correlation matrix between `bact_density_*` and `*_weighted_density` columns per strain per time point. A scatter plot of `bact_density_perinuclear` vs `inner_weighted_density` for each bacterial strain would be the first visualization to look at — do bacteria exclude actin, colocalize with it, or show no relationship?

A planned addition is a **colocalization index**: a normalized score per cell that captures how much the bacterial and actin distributions overlap across the three zones. This would be a single per-cell number easier to compare across conditions than six separate density values.

The 1h, 2h infection time points in the current dataset are the ideal testbed — we expect to see a shift from mostly cortical bacteria (early invasion, at the membrane) to more perinuclear bacteria (later, intracellular trafficking toward the perinuclear space) as infection progresses.

### 5. Image Data Organization and BRC NAS Migration

The current data is spread across local machines and scratch directories (`/scratch/sjuhasz/ActinTopologic/`). Before scaling up batch processing, the raw images need to be organized for long-term storage on the BRC (Biological Research Centre) centralized NAS system.

Proposed directory structure for the NAS:
```
BRC_NAS/
└── ActinTopology/
    └── {CellLine}_{Strain}_{Date}/
        ├── raw/
        │   └── {PlateID}/
        │       └── {WellID}/
        │           └── {FieldID}.tif
        ├── renamed_images/
        │   └── Z_Stack_Merged/
        │       └── *.tif   (pipeline-ready, 4-channel merged)
        └── results/
            ├── cytoskeleton_analysis_v3_{name}_{timestamp}.csv
            ├── {name}_plots.pdf
            └── config_used.json  (exact config snapshot used for this run)
```

Key decisions to make before migration:
- Naming convention needs to be locked in (currently inconsistent between experiments)
- Whether to keep raw Z-slices or only merged Z-stack TIFFs on NAS (raw is ~10× larger)
- Metadata file per experiment (microscope settings, objective, pixel size) should accompany each dataset

The config snapshot (`config_used.json`) is important for reproducibility — currently the same config template gets reused and modified, making it hard to know exactly what settings were used for a past run. The pipeline already saves a timestamp in the CSV filename; adding a config copy to the results directory would complete the audit trail.

### 6. PULI Cluster Integration at HCEMM

The local test runs on single experiments are fast enough (minutes per field), but a full batch across 4+ experiments with hundreds of fields each will take hours on a single GPU workstation. The PULI cluster at HCEMM has GPU nodes that would bring this down significantly.

The pipeline is already designed for this — `main_pipeline_v3.py` reads a config file and processes everything sequentially, which is easy to submit as a single cluster job. The main things to sort out:

1. **Environment setup**: The `cellpose_gpu` conda environment needs to be replicated on PULI nodes (or a container built). All dependencies (cellpose, scikit-image, tifffile, scipy, pandas) are standard.

2. **Data access**: PULI jobs need access to the NAS data — this means either mounting the NAS on compute nodes or staging data to local scratch before the job runs. The output CSVs are small and can be transferred back after the run.

3. **Job script**: A SLURM job script that activates the environment, sets `PYTHONIOENCODING=utf-8`, and calls `python main_pipeline_v3.py config_batch.json`. Estimated runtime on PULI GPU nodes: 2–4 hours for a full 4-experiment batch vs. 12–16 hours on local workstation.

4. **Parallelism**: Currently experiments run sequentially within the batch. A future improvement is to submit each experiment as a separate SLURM array job, so 4 experiments run in parallel. This requires splitting the batch config into per-experiment configs, which the pipeline already supports (just run `main_pipeline_v3.py config_experiment_A.json` four times in parallel).

---

*Written February 2026 | TopoActin V3.3D development log*
