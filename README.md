# Actin ring measurement: TopoActin V3 analysis guide (BasR vs KpSen, A549, 1 h)

This repository holds the **TopoActin V3** cytoskeleton pipeline as it stood in January–March 2026,
and the scripts and notebook needed to reproduce the paper's actin measurement:

> Compared with the ancestor strain (KpSen), basR G53V raises the **Sato-weighted skeleton density of the
> cortical actin** (`cortical_weighted_density`, outer third of the cell) in A549 cells after 1 h of infection.

**The measurement.** For every pixel of a cell, r = d_nucleus / (d_nucleus + d_membrane) runs from 0 at the
nucleus to 1 at the membrane. Actin fibres are enhanced with the multiscale Sato tubeness filter, thresholded
and thinned to a one-pixel skeleton. In each radial bin (or zone), the Sato response is summed over the
skeleton pixels and divided by the bin's area, so thick, bright bundles count more than faint filaments
and bins of different size stay comparable. Zones: inner r < 0.33, mid 0.33–0.66, outer (cortical) ≥ 0.66.
`--rings` and the notebook also give 10 rings (ring 0 = nucleus, rings 1–9 = cytoplasm) and the rings
merged 3-3-3.

The guide covers four things, from quickest to slowest:

1. **Statistics only.** Re-do the paper's statistics from a per-cell result table (minutes, no images).
2. **One field, step by step.** `notebooks/single_image_walkthrough.ipynb` takes one field (by default
   C9_F13, the BasR field of the paper figure) from the raw stack to the paper figure (a few minutes).
3. **Representative images.** Re-run the image analysis on the 18 fields the paper's example images
   come from, and compare cell by cell with the original results (about 30 min on a GPU).
4. **Whole plates.** Re-run a full plate (about 600 fields) and repeat the statistics on it (hours).

**Example data to test with:** `examples/` holds four fields of the paper plate Exp10 (BasR C9_F13 and
D10_F02, KpSen B5_F11, control D11_F09) and their January 2026 result rows (see `examples/README.md`).
Whole plates are not stored here; for Steps 3–4 you need the plate folders (see
[What you need](#what-you-need)).

---

## Contents

| Path | What it is |
|---|---|
| `old_pipeline/` | The V3 code, unchanged. `main_pipeline_v3.py` is the batch pipeline; `TopoActin_V3_Interactive_DOWNSAMPLING.ipynb` is the interactive notebook for looking at one image step by step. `docs/` has the original design notes (`ARCHITECTURE.md`, `WORKBOOK_2026-02.md`). |
| `examples/` | Four example fields (merged Z-stacks) and their January 2026 result rows, for testing. |
| `notebooks/single_image_walkthrough.ipynb` | One field from the raw stack to the paper figure: segmentation, Sato filter, skeleton, radial zones, 10 rings, figure. |
| `replication/v3_paper.py` | Applies the paper settings to the unchanged V3 code at run time; runs one field with the V3 classes; the 10-ring profile. |
| `replication/notebook_method.py` | The per-cell measurement of the interactive V3 notebook (Steps 6b/6c), which made the paper figure. |
| `replication/ring_plots.py` | Figures for the walkthrough notebook, in the style of the paper figure. |
| `replication/config_paper.json` | The settings to use (see [Settings](#settings-used)). |
| `replication/field_lists/` | Which fields to analyse: the exact fields of each paper plate, the 18 representative fields, and the held-out check fields. |
| `replication/build_stacks.py` | Rebuilds `Z_Stack_Merged/` from a raw Harmony export, byte-identical to the stacks the paper used. |
| `replication/prepare_input.py` | Puts the chosen fields of a plate into one input folder. |
| `replication/run_v3.py` | Runs the old pipeline on that folder, with intensities measured the way the paper runs did. |
| `replication/paper_stats.py` | The paper's statistics (fidelity filter, winsorizing, well medians, Mann-Whitney, bootstrap, BH). |
| `replication/compare_to_reference.py` | Compares a re-run with the original result table. |
| `replication/plot_figure_FG.py` | Redraws paper Fig. F–G (well-median violins) from `paper_stats.py` output. |
| `environment.yml` | Python environment. |
| `LICENSE` | MIT. |

---

## What you need

**Images.** One folder per plate, as exported from Harmony and merged into Z-stacks. The scripts
use the `Z_Stack_Merged/` subfolder, one multi-channel TIFF per field, named
`<Well>_F<nn>_Z_stack_merged.tiff` (e.g. `C9_F13_Z_stack_merged.tiff`).
Channel order: **actin, nuclei, brightfield, bacteria** (the paper plates have no bacteria stain).

The three paper plates:

| Paper ID | Plate folder | Fields per well on disk | Fields analysed |
|---|---|---|---|
| Exp1 | `2025-11-11_A549_Cellmask_KpsenBasR_1h_plate2_2025-11-12T10_43_01-Measurement 1` | 26 | 10 per well, 599 in total |
| Exp8 | `2026-01-13_A549_Cellmask_KpSenBasR_1h_plate2_2026-01-14T13_24_13-Measurement 1` | 14 | 10 per well, 597 in total |
| Exp10 | `2026-01-13_A549_Cellmask_KpSenBasR_1h_plate4_2026-01-14T12_48_47-Measurement 1` | 14 | 10 per well, 597 in total |

Each plate folder also has a `results/` subfolder holding the **original January 2026 result
table** (`…Measurement 1.csv`, one row per cell). That table is the reference you compare with.
All of the paper's representative images come from **Exp10**.

**No `Z_Stack_Merged/` folder?** If you only have the raw Harmony export (`Images/` with
`Index.xml` and one TIFF per plane), rebuild the stacks first:

```bash
python replication/build_stacks.py "<plate folder>"
```

The original stacks were made in ImageJ 1.54p, and they hold the raw planes unchanged (checked plane by
plane). `build_stacks.py` picks the channels by name from `Index.xml` and puts actin (CellMask)
first and nuclei (Hoechst) second, the order V3 expects. It reproduced all 18 representative stacks
byte for byte. Do not use a converter that orders channels by Harmony channel number: on these plates
that puts Hoechst first, and V3 would read the nuclei as actin.

**Plate layout** (same on every plate; groups are set by the column number):

| Columns | Group |
|---|---|
| 2 and 11 | Control (uninfected) |
| 3, 4, 5, 6 | KpSen (ancestor) |
| 7, 8, 9, 10 | BasR (basR G53V) |

Rows B–G, so 24 KpSen wells, 24 BasR wells and 12 Control wells. If your plate uses a different
layout, edit `well_to_group` in `replication/config_paper.json`.

**Computer.** An NVIDIA GPU is strongly recommended. One field takes about 1–1.5 min on a small
4 GB card. Without a GPU, add `--cpu` (expect it to be many times slower).

---

## Step 0: install

```bash
git clone https://github.com/Sompaa/actin-ring-measurement.git
cd actin-ring-measurement
conda env create -f environment.yml
conda activate topoactin_v3
```

If you already have a Cellpose 4 environment (e.g. `topoactin_sam` on the lab machines), you can use
that instead. Cellpose downloads its `cpsam` model the first time it runs (about 1 GB).

All commands below are run from the repository folder.

---

## Step 1: statistics from an existing result table

Start here, because it needs no images. It answers "does this table show the paper result?"

```bash
python replication/paper_stats.py \
    "<Exp1 plate>/results/<…>.csv" "<Exp8 plate>/results/<…>.csv" "<Exp10 plate>/results/<…>.csv" \
    --names Exp1 Exp8 Exp10 -o stats_paper/
```

What it does, for each plate separately (this is the method written on the paper's statistics slide):

1. Keep cells with **fidelity ≥ 0.25** (how well the skeleton follows the actin; about 96–97% of cells pass).
2. Inside each well, pull the most extreme cells in to that well's 5th and 95th percentile
   (winsorizing). Outlier cells are kept but cannot drag the well's value around.
3. Take **one median per well** for every measurement. The well is the real replicate; cells in the
   same well are not independent of each other.
4. Compare the 24 KpSen wells with the 24 BasR wells with a **Mann-Whitney U test** (a rank test
   that does not assume a bell-shaped distribution).
5. Report the **difference of the medians** (BasR minus KpSen) in the measurement's own units, with
   a 95% range from 2000 bootstrap repeats (the wells are re-drawn at random 2000 times to see how
   much the difference could wobble).
6. Correct for having tested ~39 measurements at once (**Benjamini-Hochberg**, column `q_BH`).

Output: `stats_<name>.csv` (every measurement) and `stats_headline.csv` (the two paper measurements).

With two or more tables it also runs the test from the final paper figure: all plates pooled, one
value per well, a **linear mixed-effects model** with the experiment as a random effect (so
plate-to-plate differences in overall brightness do not count as a treatment effect), then
Benjamini-Hochberg. That writes `stats_pooled_lmm.csv` and `well_medians_pooled.csv`. Redraw the
figure panels with:

```bash
python replication/plot_figure_FG.py stats_paper/
```

**What you should see on the original tables** (checked 2026-10-08):

| Plate | Measurement | KpSen | BasR | BasR − KpSen (95% range) | q (BH) |
|---|---|---|---|---|---|
| Exp1 | cortical_intensity | 12,290 | 15,550 | +3,269 (2,135 to 4,154) | 0.0001 |
| Exp1 | cortical_weighted_density | 164.2 | 191.3 | +27.1 (6.1 to 47.5) | 0.021 |
| Exp8 | cortical_intensity | 10,470 | 12,970 | +2,495 (585 to 3,624) | 0.007 |
| Exp8 | cortical_weighted_density | 148.5 | 173.7 | +25.1 (10.5 to 50.2) | 0.007 |
| Exp10 | cortical_intensity | 10,180 | 13,590 | +3,409 (1,639 to 6,270) | 0.002 |
| Exp10 | cortical_weighted_density | 150.0 | 196.6 | +46.6 (22.2 to 68.6) | 0.003 |

BasR is higher on both measurements in all three experiments, and every difference stays
significant after correction. So the statistics step reproduces the paper claim from the original
tables. (The bootstrap range can move slightly between computers; the medians and p-values do not.)

Pooled mixed model over the three plates: `cortical_intensity` +3,285 (95 % CI 2,538 to 4,031,
q = 2×10⁻¹⁶), `cortical_weighted_density` +33.8 (23.7 to 44.0, q = 1×10⁻⁹), both "****" as in the
paper. The well-median ranges equal Fig. F–G exactly: ancestor 5,817–15,763 / BasR 8,811–21,117
(intensity), 90–230 / 120–277 (bundle density); control medians 10,911 and 145.

One thing to keep in mind when reading the full tables: the number of measurements that differ
varies a lot between plates (7 of 39 in Exp1, 14 in Exp8, 32 in Exp10), so "only the cortex
changes" holds best on Exp1.

---

## Step 2: one field, step by step (the paper figure)

`notebooks/single_image_walkthrough.ipynb` analyses one field and ends with the paper figure. By default
it uses the example `examples/C9_F13_Z_stack_merged.tiff`; run it top to bottom:

```bash
jupyter lab notebooks/single_image_walkthrough.ipynb
```

For another field set `IMAGE` in the first cell (or `ACTIN_IMAGE=<stack> jupyter lab ...`). It shows, in order: the Cellpose segmentation and the cells kept by
the filters; the Sato tubeness filter on the field; the measurement of every kept cell; for one chosen
cell the eight intermediate panels (cell, Sato, threshold, skeleton, zones, radial profile, numbers); the
three panels of the paper figure (original image, radial overlay, distribution plot); and the 10-ring
profile with the rings merged 3-3-3. Figures (PNG at 300 dpi, SVG) and the per-cell table go to
`walkthrough_out/`.

**Check:** on the example C9_F13, cell 21 gives inner / mid / outer Sato-weighted density
0.0015 / 0.0028 / 0.0048, the values printed for that cell on the February 2026 figure slide, with
the same skeleton and radial profile (tested 2026-10-13, Cellpose 4.2.1.1, GTX 1650; about 1 min).

**Which method makes the figure.** The paper figure was drawn with the interactive V3 notebook
(Steps 6b/6c of `old_pipeline/TopoActin_V3_Interactive_DOWNSAMPLING.ipynb`), not with the batch run.
`replication/notebook_method.py` is that per-cell method taken out of the notebook's widgets, unchanged:
Sato (scikit-image defaults) on the 0–1 rescaled actin maximum projection, threshold at the 55th
percentile of the cell's Sato values, objects < 3 px removed, closing radius 1, skeleton, r from the
cell's own nucleus, 10 bins of r. The densities are shown unscaled, as on the paper figure
(0.001–0.007). Later versions of the notebook (from mid-February 2026) multiplied the displayed
density by 5 (and intensities by 10) only to put the numbers on the same scale as the batch tables;
the multiplier is not part of the measurement (`multiplier=5` reproduces it).

The statistics (Step 1) use the batch tables (`main_pipeline_v3.py`); the figure shows the same
quantity for one cell, measured with the notebook's per-cell settings.

---

## Step 3: the representative images

The paper's example images come from 25 fields of Exp10. 18 of them were also in the analysed set,
so they have original numbers to compare with; they are listed in
`replication/field_lists/Exp10_representative_fields.txt`. (The other 7, e.g. D10_F02 and B5_F14,
were shown in figures but never measured: `Exp10_representative_not_analysed.txt`.)

**2a. Collect the fields** (makes links, nothing is copied; add `--copy` on Windows):

```bash
python replication/prepare_input.py "<Exp10 plate folder>" \
    replication/field_lists/Exp10_representative_fields.txt  input_rep/
```

**2b. Run the pipeline:**

```bash
python replication/run_v3.py input_rep/ out_rep/
```

The per-cell table appears in `out_rep/results/cytoskeleton_analysis_v3_Single_Experiment_<date>.csv`,
with skeleton images and comparison figures for the first cells beside it. `out_rep/settings_used.json`
records exactly what was run.

**2c. Compare with the original:**

```bash
python replication/compare_to_reference.py out_rep/results/cytoskeleton_analysis_v3_*.csv \
    "<Exp10 plate folder>/results/<…>.csv" -o compare_rep/
```

You get `cells_per_field.csv` (how many cells each run kept), `metric_agreement.csv` (for cells
found by both runs: the typical ratio new/original and how closely they track each other) and
`well_medians.csv`. Cells are paired by size and shape, because a cell's number (`cell_id`) is not
the same between two runs.

**To look at one image interactively**, open `old_pipeline/TopoActin_V3_Interactive_DOWNSAMPLING.ipynb`
in Jupyter, change `img_dir` in the "STEP 1: LOAD IMAGE" cell to your `input_rep/` folder (it points
at an old Windows path), and run the cells top to bottom. Each step has sliders and a preview: brightness
and contrast, Cellpose segmentation, nucleus and cell filters, Sato filter + skeleton, the three
radial zones, and a 4-panel figure. The notebook writes no result table; it is for looking and for
choosing settings.

See [What to expect](#what-to-expect-from-a-re-run) for how close a re-run gets today.

---

## Step 4: whole plates

Same three commands, with the full field list of the plate instead of the 18 fields:

```bash
python replication/prepare_input.py "<Exp10 plate folder>" replication/field_lists/Exp10_paper_fields.txt input_exp10/
python replication/run_v3.py input_exp10/ out_exp10/
python replication/compare_to_reference.py out_exp10/results/cytoskeleton_analysis_v3_*.csv \
    "<Exp10 plate folder>/results/<…>.csv" -o compare_exp10/
python replication/paper_stats.py out_exp10/results/cytoskeleton_analysis_v3_*.csv --names Exp10_rerun -o stats_exp10/
```

Repeat with `Exp1_paper_fields.txt` / `Exp8_paper_fields.txt` for the other plates. Use `all`
instead of a list file to analyse every field on the plate (Exp1 has 26 per well, so ~1560 fields).

About 600 fields take roughly 10–15 h on a 4 GB GPU; run it overnight or on the cluster. On a whole
plate, `well_medians.csv` is the most useful comparison: plot new against original per well, and the
points should lie close to the diagonal.

**Why the field lists, and not the pipeline's own downsampling?** The paper runs analysed 10 random
fields per well (the same 10 in every well). The `position_downsampling.py` module with its seed 42
does **not** pick those same fields again on these plates (about half differ), probably because the
original run listed files in a different order on Windows. The lists in `field_lists/` are read from
the original result tables, so they are exactly the paper's fields. Keep `downsampling.enabled` off.

---

## Settings used

The settings that produced the January 2026 tables were never saved, so they were reconstructed.
The rule was strict: **only values with written evidence were allowed** (a hard edge in the January
tables, a saved config, a value printed in a notebook output, or a statement from the people who ran
it). Among those, the choice was made on the 18 Exp10 representative fields, then checked on 20 fields
of Exp1 and Exp8 that were picked at random *before* calibration and never used to choose anything
(`field_lists/Exp*_holdout_fields.txt`).

`replication/config_paper.json` holds the result. Compared with the February batch config
(`old_pipeline/config_batch.json`), these values changed:

| Setting | Feb config | Used here | Evidence |
|---|---|---|---|
| 3D mode | off | on (`keep_3d_volume`) | the paper tables say `processing_mode = 3D` |
| Actin brightness / contrast | 1.4 / 1.2 | **1.4 / 1.75** | worked back from January's own `cortical_intensity` and `inner_intensity` (the adjusted image is a fixed formula of the raw image, so the setting can be solved from the measured values; best-fitting region brightness × contrast ≈ 2.45). 1.4 is the batch-config brightness, 1.75 the contrast slider printed in the 17 Jan notebook |
| Intensities measured on | raw image (since 10 Feb) | adjusted image | raw gives values 12–25× too low; the January runs predate the 10 Feb change |
| Cell max eccentricity | 0.85 | **0.95** | hard edge in the January tables: max = 0.950 over 30,657 cells, 17 % above 0.85 |
| Cell min solidity | 0.70 | 0.70 | hard edge: min = 0.700 |
| Cell min area | 3000 px | **500 px** | January kept cells down to 647 px (1 % below 3000); 500 is the V2 value |
| Border buffer | 20 px | 20 px | loosening it keeps 127 extra border cells that January never kept |
| Nucleus max eccentricity | 0.80 | **0.85** | V2 (14 Jan) config; 0.80 alone drops 31 of January's cells |
| Nucleus proximity exclusion | on | **off** | not present in the January-era code |
| Sato filter | scikit-image default (dark ridges) | **bright ridges** (`sato_black_ridges: false`) | the dark-ridge default traces the gaps between fibres; the predecessor tool (Skeletactin, Oct 2025) used bright ridges |
| Sato σ | 0.5–2.0, 5 steps | **0.7, 0.96, 1.22, 1.48, 1.74, 2.0** | January values (V2 config, code default, confirmed by the authors) |
| Min object size / closing | 60 px / 0 | **3 px / 1** | slider values printed in the 17 Jan notebook and in the 17 Feb "Old pipeline" reference block |
| Extra Vieira smoothing | applied (bug: σ 0.5 ran although disabled) | off | disabled in every January-era setting |

Unchanged: nuclei brightness/contrast 1.4 / 1.3, Cellpose `cpsam` with automatic diameter, flow 0.4,
cell probability 0.0; Gaussian σ 0.7 before Sato; Li threshold + hysteresis (low 0.4); zones inner
0–0.33, mid 0.33–0.66, cortical 0.66–1.0 of the nucleus-to-edge distance; cortex = outer 5 px.

Tested and rejected, with the reason: the 17 Jan notebook's Cellpose sliders (flow 0.25, probability
−1) and its brightness/contrast 1.9 / 1.75 (fewer of January's cells, intensities ×1.3–1.5 too high);
the custom 3D Sato from the notebooks (no cell-by-cell agreement); a 55th-percentile threshold (kept as
a sensitivity check: weighted densities ×1.12 instead of ×1.21, but branch counts ×0.72).

**Cellpose model.** The V3 code asks for `pretrained_model='cpsam'` (Cellpose-SAM, Cellpose 4), the
pretrained model used as published, without retraining. The exact Cellpose version of the January runs
was not recorded, and the lab's Windows environment was later downgraded to Cellpose 3.1, where the name
`cpsam` does not exist and Cellpose silently falls back to `cyto3`. So it was checked on the 18
representative Exp10 fields (same settings, only the model swapped):

| Model | Cells kept | Median cell area (px) | January cells found again |
|---|---|---|---|
| January original | 252 | 19,704 | – |
| Cellpose 3, as coded (falls back to `cyto3`) | 303 | 6,153 | 0 |
| `cyto2` cells + `nuclei` | 97 | 9,330 | 6 |
| `cyto2` cells + `cyto2` nuclei | 118 | 9,128 | 5 |
| **Cellpose 4 `cpsam`** | **256** | **20,170** | **211 (84 %)**, cell area r = 0.999 |

Only `cpsam` gives January's cells, so the paper runs used Cellpose 4 `cpsam`.

`run_v3.py` applies the two code-level switches (intensity on the adjusted image, bright-ridge Sato)
at run time through `replication/v3_paper.py`; the files in `old_pipeline/` are not edited. Every run writes `settings_used.json`.

---

## What to expect from a re-run

Measured 2026-10-09 with `config_paper.json` through `run_v3.py`, compared with the original January
tables by `compare_to_reference.py`. "Found again" = January cells that the re-run also keeps (paired by
size and shape). Ratios are medians of new ÷ January over those cells; r is the cell-by-cell correlation.

| | Exp10, 18 fields (calibration) | Exp1, 10 fields (held out) | Exp8, 10 fields (held out) |
|---|---|---|---|
| Cells kept: January / re-run | 252 / 256 | 119 / 120 | 218 / 209 |
| January cells found again | 211 (84 %) | 110 (92 %) | 167 (77 %) |
| cell_area | 1.00 (r 0.999) | 1.00 (r 0.999) | 1.01 (r 0.999) |
| **cortical_intensity** | **1.07 (r 0.96)** | **1.11 (r 0.99)** | **1.06 (r 0.88)** |
| inner_intensity | 1.03 (r 0.96) | 1.08 (r 0.99) | 1.08 (r 0.92) |
| **cortical_weighted_density** | **1.20 (r 0.84)** | **0.98 (r 0.89)** | **1.15 (r 0.82)** |
| skeleton_density | 1.11 | 1.11 | 1.16 |
| fidelity | 1.30 | 1.30 | 1.39 |

Before calibration (February settings) the same Exp10 fields gave 183 cells with 56 % found again,
`cortical_intensity` ×0.88 and `cortical_weighted_density` ×0.70, with branch counts ×1.5.

What this means:

- **The same cells.** 77–92 % of January's cells come out again, and about as many cells are kept
  in total. Most of the rest is Cellpose: for ~8 % of January's cells there is no matching mask at all
  today (the exact Cellpose 4 release of January was not recorded; see [Settings](#settings-used)).
- **The paper's two measurements are on the same scale and track January cell by cell**, within
  2–20 % depending on the plate.
- **Fidelity reads higher** (×1.3–1.4) and does not track January cell by cell; every re-run cell
  passes the ≥ 0.25 quality cut, against ~97 % in January. The skeleton's fine structure is the part
  that could not be fully reconstructed. Branch counts do match (×1.04, r 0.94, on Exp10).
- **Per well**, with only 2 fields per held-out well, most wells are within ~10 % of January and a
  few differ by 20–40 %. Whole plates (10 fields per well) are the real test: run Step 4 and plot
  `well_medians.csv` new against January.

For the paper claim, what must hold on the re-run plates is the **direction and significance** of
BasR − KpSen in Step 1's statistics (both the per-plate test and the pooled mixed model), and well
medians in the same range as Fig. F–G.

---

## Troubleshooting

- **`No *_Z_stack_merged.tif(f) files found`**: point `prepare_input.py` at the plate folder or at its
  `Z_Stack_Merged/` subfolder. Raw Harmony `Images/` are single planes and are not read directly.
- **`WARNING: … listed field(s) are not on this plate`**: you used the field list of another plate.
- **Cellpose falls back to the CPU** although you have a GPU: your PyTorch build has no CUDA. Reinstall
  `torch` from pytorch.org with the CUDA version of your driver.
- **Out of GPU memory** on a small card: close other GPU programs; one field at a time needs about 3.4 GB.
- **Intensities 10–25× lower than the original**: you ran with `--measure-on raw`, or called
  `old_pipeline/main_pipeline_v3.py` directly. Use `replication/run_v3.py`.
- **All cells end up in group "Unassigned"**: the file names do not start with `<Row><Column>_`
  (e.g. `C9_`), or the plate layout differs; edit `well_to_group`.
