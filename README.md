# Actin ring measurement

Code, example data and instructions for the actin cytoskeleton measurement in A549 cells infected for
1 h with *Klebsiella pneumoniae* KpSen (ancestor) or the basR G53V mutant (BasR). The analysis is the
TopoActin V3 pipeline in the state used for the original analysis (January to March 2026).

The reported result is that BasR increases the Sato-weighted skeleton density of cortical actin
(`cortical_weighted_density`, outer third of the cell) relative to KpSen.

## Method

For each pixel of a segmented cell, the normalized radial position
r = d_nucleus / (d_nucleus + d_membrane) is 0 at the nuclear boundary and 1 at the plasma membrane.
F-actin is enhanced with the multiscale Sato tubeness filter, thresholded and reduced to a one-pixel
skeleton. For each radial bin or zone, the Sato response summed over the skeleton pixels is divided by
the area of the bin. Thick, bright bundles therefore contribute more than faint filaments, and bins of
different area can be compared directly. Zones are defined as inner (r < 0.33), mid (0.33 ≤ r < 0.66)
and outer, or cortical (r ≥ 0.66). The walkthrough notebook and `run_v3.py --rings` also report ten
rings (ring 0 is the nucleus, rings 1 to 9 divide the cytoplasm into equal steps of r) and the rings
merged in groups of three.

Cells and nuclei are segmented with the pretrained Cellpose-SAM model (`cpsam`, Cellpose 4), used
without retraining.

## Repository contents

| Path | Content |
|---|---|
| `old_pipeline/` | TopoActin V3 source code, unmodified. `main_pipeline_v3.py` is the batch pipeline; `TopoActin_V3_Interactive_DOWNSAMPLING.ipynb` is the interactive notebook used for single-image inspection and for the figure. `docs/` contains the original design notes (`ARCHITECTURE.md`, `WORKBOOK_2026-02.md`). |
| `examples/` | Four example fields (merged Z-stacks) with their original per-cell results. |
| `notebooks/single_image_walkthrough.ipynb` | Analysis of one field from the raw stack to the figure: segmentation, Sato filter, skeleton, radial zones, ten rings. |
| `replication/v3_paper.py` | Applies the analysis settings to the unmodified V3 code at run time, analyses a single field with the V3 classes, and computes the ten-ring profile. |
| `replication/notebook_method.py` | Per-cell measurement of the interactive V3 notebook (Steps 6b and 6c), which produced the figure. |
| `replication/ring_plots.py` | Plotting functions for the walkthrough notebook. |
| `replication/config_paper.json` | Analysis settings (see [Analysis settings](#analysis-settings)). |
| `replication/field_lists/` | Fields analysed on each plate, the 18 representative fields, and the held-out validation fields. |
| `replication/build_stacks.py` | Builds `Z_Stack_Merged/` from a Harmony export; output is byte-identical to the stacks used in the original analysis. |
| `replication/prepare_input.py` | Collects selected fields of a plate into one input folder. |
| `replication/run_v3.py` | Runs the V3 batch pipeline on an input folder with the analysis settings. |
| `replication/paper_stats.py` | Statistical analysis (fidelity filter, winsorization, well medians, Mann-Whitney U test, bootstrap confidence intervals, Benjamini-Hochberg correction, linear mixed model). |
| `replication/compare_to_reference.py` | Cell-by-cell and well-level comparison of a re-run with the original results. |
| `replication/plot_figure_FG.py` | Well-median violin plots (figure panels F and G) from the output of `paper_stats.py`. |
| `environment.yml` | Conda environment. |
| `LICENSE` | MIT licence. |

## Data

The example fields in `examples/` are sufficient to run the walkthrough notebook and to test the batch
pipeline (see `examples/README.md`). Re-analysis of representative fields or whole plates requires the
plate folders, which are not included in this repository.

Each plate folder contains a `Z_Stack_Merged/` subfolder with one multichannel TIFF per field, named
`<Well>_F<nn>_Z_stack_merged.tiff` (for example `C9_F13_Z_stack_merged.tiff`). Channel order is actin
(CellMask), nuclei (Hoechst), brightfield, bacteria; the plates analysed here have no bacterial stain.
The `results/` subfolder of each plate holds the original per-cell result table (`…Measurement 1.csv`),
which serves as the reference for comparisons.

| Experiment | Plate folder | Fields per well acquired | Fields analysed |
|---|---|---|---|
| Exp1 | `2025-11-11_A549_Cellmask_KpsenBasR_1h_plate2_2025-11-12T10_43_01-Measurement 1` | 26 | 10 per well, 599 in total |
| Exp8 | `2026-01-13_A549_Cellmask_KpSenBasR_1h_plate2_2026-01-14T13_24_13-Measurement 1` | 14 | 10 per well, 597 in total |
| Exp10 | `2026-01-13_A549_Cellmask_KpSenBasR_1h_plate4_2026-01-14T12_48_47-Measurement 1` | 14 | 10 per well, 597 in total |

All representative images were acquired on Exp10.

If only the raw Harmony export is available (`Images/` with `Index.xml` and one TIFF per plane), build
the stacks first:

```bash
python replication/build_stacks.py "<plate folder>"
```

The original stacks were assembled in ImageJ 1.54p from the unmodified raw planes. `build_stacks.py`
selects channels by name from `Index.xml` and writes actin (CellMask) first and nuclei (Hoechst)
second, the order V3 expects; it reproduces all 18 representative stacks byte for byte. Converters
that order channels by Harmony channel number place Hoechst first and must not be used.

Plate layout (identical on all plates; rows B to G):

| Columns | Group |
|---|---|
| 2 and 11 | Control (uninfected) |
| 3, 4, 5, 6 | KpSen (ancestor) |
| 7, 8, 9, 10 | BasR (basR G53V) |

This gives 24 KpSen, 24 BasR and 12 control wells per plate. For a different layout, edit
`well_to_group` in `replication/config_paper.json`.

Hardware: an NVIDIA GPU is recommended. One field takes 1 to 1.5 min on a 4 GB GPU. Add `--cpu` to run
without a GPU, at a considerably lower speed.

## Installation

```bash
git clone https://github.com/Sompaa/actin-ring-measurement.git
cd actin-ring-measurement
conda env create -f environment.yml
conda activate topoactin_v3
```

Any environment with Cellpose 4 can be used instead. The `cpsam` model weights (about 1 GB) are
downloaded on first use. All commands below are run from the repository root.

## 1. Statistical analysis from per-cell tables

This step requires no images.

```bash
python replication/paper_stats.py \
    "<Exp1 plate>/results/<…>.csv" "<Exp8 plate>/results/<…>.csv" "<Exp10 plate>/results/<…>.csv" \
    --names Exp1 Exp8 Exp10 -o stats_paper/
```

For each plate separately the script:

1. keeps cells with skeleton fidelity ≥ 0.25 (96 to 97 % of cells);
2. winsorizes each measurement within each well at the 5th and 95th percentiles;
3. takes the median of each measurement per well, the well being the experimental unit;
4. compares the 24 KpSen and 24 BasR wells with a two-sided Mann-Whitney U test;
5. reports the difference of medians (BasR minus KpSen) with a 95 % confidence interval from 2,000
   bootstrap resamples of wells;
6. corrects for the approximately 39 measurements tested with the Benjamini-Hochberg procedure
   (column `q_BH`).

Output: `stats_<name>.csv` (all measurements) and `stats_headline.csv` (the reported measurements).

With two or more tables, the script also fits a linear mixed-effects model to the pooled well medians
with experiment as a random effect, followed by Benjamini-Hochberg correction. Between-plate
differences in overall intensity are thereby not attributed to the infection. Results are written to
`stats_pooled_lmm.csv` and `well_medians_pooled.csv`. Figure panels F and G are drawn with:

```bash
python replication/plot_figure_FG.py stats_paper/
```

Results on the original tables (computed 2026-10-08):

| Plate | Measurement | KpSen | BasR | BasR − KpSen (95 % CI) | q (BH) |
|---|---|---|---|---|---|
| Exp1 | cortical_intensity | 12,290 | 15,550 | +3,269 (2,135 to 4,154) | 0.0001 |
| Exp1 | cortical_weighted_density | 164.2 | 191.3 | +27.1 (6.1 to 47.5) | 0.021 |
| Exp8 | cortical_intensity | 10,470 | 12,970 | +2,495 (585 to 3,624) | 0.007 |
| Exp8 | cortical_weighted_density | 148.5 | 173.7 | +25.1 (10.5 to 50.2) | 0.007 |
| Exp10 | cortical_intensity | 10,180 | 13,590 | +3,409 (1,639 to 6,270) | 0.002 |
| Exp10 | cortical_weighted_density | 150.0 | 196.6 | +46.6 (22.2 to 68.6) | 0.003 |

BasR exceeds KpSen for both measurements in all three experiments, and all differences remain
significant after correction. Bootstrap intervals can vary slightly between machines; medians and
p-values are deterministic.

Pooled mixed model over the three plates: `cortical_intensity` +3,285 (95 % CI 2,538 to 4,031,
q = 2×10⁻¹⁶), `cortical_weighted_density` +33.8 (23.7 to 44.0, q = 1×10⁻⁹). The well-median ranges
match figure panels F and G: KpSen 5,817 to 15,763 and BasR 8,811 to 21,117 (intensity), 90 to 230 and
120 to 277 (weighted density); control medians 10,911 and 145.

The number of measurements that differ between KpSen and BasR varies between plates (7 of 39 in Exp1,
14 in Exp8, 32 in Exp10). A change restricted to the cortex is therefore best supported by Exp1.

## 2. Single-field analysis and figure

`notebooks/single_image_walkthrough.ipynb` analyses one field and produces the figure. By default it
uses `examples/C9_F13_Z_stack_merged.tiff`:

```bash
jupyter lab notebooks/single_image_walkthrough.ipynb
```

To analyse another field, set `IMAGE` in the first cell or start Jupyter with `ACTIN_IMAGE=<stack>`.
The notebook shows the Cellpose segmentation and the cells passing the filters, the Sato response of
the field, the measurements of all retained cells, the intermediate steps for one selected cell (cell
region, Sato response, threshold, skeleton, zones, radial profile), the figure panels (original image,
radial overlay, distribution plot) and the ten-ring profile. Figures (PNG at 300 dpi and SVG) and the
per-cell table are written to `walkthrough_out/`.

Validation: for cell 21 of the example C9_F13, the notebook gives inner, mid and outer Sato-weighted
densities of 0.0015, 0.0028 and 0.0048, identical to the values of the original analysis for this
cell, with the same skeleton and radial profile (Cellpose 4.2.1.1, GTX 1650, about 1 min per field).

The figure was produced with the interactive V3 notebook (Steps 6b and 6c of
`old_pipeline/TopoActin_V3_Interactive_DOWNSAMPLING.ipynb`), not with the batch pipeline.
`replication/notebook_method.py` contains this per-cell procedure without the notebook widgets and
without changes: Sato filter (scikit-image defaults) on the actin maximum projection rescaled to 0 to 1,
threshold at the 55th percentile of the Sato response within the cell, removal of objects smaller than
3 px, morphological closing with radius 1, skeletonization, radial position relative to the cell's own
nucleus, and ten bins of r. Densities are unscaled (0.001 to 0.007), as in the figure. Later versions
of the notebook multiplied the displayed density by 5 and intensities by 10 to match the value range
of the batch tables; this factor is a display scaling and not part of the measurement
(`multiplier=5` reproduces it).

The statistical analysis (section 1) uses the per-cell tables of the batch pipeline
(`main_pipeline_v3.py`). The figure shows the same quantity for individual cells, measured with the
per-cell settings of the notebook.

## 3. Representative fields

The representative images were taken from 25 fields of Exp10. Of these, 18 were part of the analysed
set and have original results for comparison; they are listed in
`replication/field_lists/Exp10_representative_fields.txt`. The remaining 7 (for example D10_F02 and
B5_F14) were used for illustration only and were not measured
(`Exp10_representative_not_analysed.txt`).

Collect the fields (creates symbolic links; add `--copy` on Windows):

```bash
python replication/prepare_input.py "<Exp10 plate folder>" \
    replication/field_lists/Exp10_representative_fields.txt  input_rep/
```

Run the pipeline:

```bash
python replication/run_v3.py input_rep/ out_rep/
```

The per-cell table is written to `out_rep/results/cytoskeleton_analysis_v3_Single_Experiment_<date>.csv`,
together with skeleton images and comparison figures for the first cells. `out_rep/settings_used.json`
records the settings of the run.

Compare with the original results:

```bash
python replication/compare_to_reference.py out_rep/results/cytoskeleton_analysis_v3_*.csv \
    "<Exp10 plate folder>/results/<…>.csv" -o compare_rep/
```

Output: `cells_per_field.csv` (cells retained per field in each run), `metric_agreement.csv` (for cells
present in both runs: median ratio re-run/original and correlation) and `well_medians.csv`. Cell
identifiers differ between runs, so cells are paired by area and perimeter.

The interactive notebook `old_pipeline/TopoActin_V3_Interactive_DOWNSAMPLING.ipynb` can be used to
inspect a single image: set `img_dir` in the "STEP 1: LOAD IMAGE" cell to the input folder and run the
cells in order. It provides previews for brightness and contrast, Cellpose segmentation, nucleus and
cell filters, the Sato filter and skeleton, the radial zones and a four-panel figure. It does not write
a result table.

## 4. Whole plates

The same commands with the complete field list of a plate:

```bash
python replication/prepare_input.py "<Exp10 plate folder>" replication/field_lists/Exp10_paper_fields.txt input_exp10/
python replication/run_v3.py input_exp10/ out_exp10/
python replication/compare_to_reference.py out_exp10/results/cytoskeleton_analysis_v3_*.csv \
    "<Exp10 plate folder>/results/<…>.csv" -o compare_exp10/
python replication/paper_stats.py out_exp10/results/cytoskeleton_analysis_v3_*.csv --names Exp10_rerun -o stats_exp10/
```

Use `Exp1_paper_fields.txt` and `Exp8_paper_fields.txt` for the other plates, or `all` to analyse every
field (Exp1 has 26 fields per well, about 1,560 in total). About 600 fields take 10 to 15 h on a 4 GB
GPU. On whole plates, `well_medians.csv` is the most informative comparison; re-run and original well
medians should lie close to the identity line.

The original analysis used 10 randomly selected fields per well, the same positions in every well.
`position_downsampling.py` with seed 42 does not reproduce this selection on these plates (about half
of the fields differ), most likely because the original run listed the files in a different order. The
lists in `field_lists/` were read from the original result tables and contain exactly the analysed
fields. Keep `downsampling.enabled` set to false.

## Analysis settings

The settings of the original batch runs were not saved and were reconstructed from documented
evidence only: boundaries visible in the original result tables, saved configuration files, values
printed in notebook outputs, and information from the people who ran the analysis. Candidate settings
were selected on the 18 representative Exp10 fields and validated on 20 fields of Exp1 and Exp8 that
were drawn at random before calibration and not used for any choice
(`field_lists/Exp*_holdout_fields.txt`).

The resulting settings are in `replication/config_paper.json`. Differences from the February batch
configuration (`old_pipeline/config_batch.json`):

| Setting | February configuration | Used here | Evidence |
|---|---|---|---|
| 3D mode | off | on (`keep_3d_volume`) | original tables report `processing_mode = 3D` |
| Actin brightness / contrast | 1.4 / 1.2 | 1.4 / 1.75 | derived from the original `cortical_intensity` and `inner_intensity` (the adjusted image is a fixed function of the raw image; best-fitting brightness × contrast ≈ 2.45); 1.4 is the batch brightness, 1.75 the contrast printed in the notebook of 17 January |
| Image used for intensities | raw image (from 10 February) | brightness/contrast-adjusted image | the raw image gives values 12 to 25 times lower; the original runs predate the 10 February change |
| Cell maximum eccentricity | 0.85 | 0.95 | maximum in the original tables is 0.950 over 30,657 cells; 17 % exceed 0.85 |
| Cell minimum solidity | 0.70 | 0.70 | minimum in the original tables is 0.700 |
| Cell minimum area | 3000 px | 500 px | original tables contain cells down to 647 px (1 % below 3000); 500 is the V2 value |
| Border exclusion | 20 px | 20 px | a smaller margin retains 127 border cells absent from the original tables |
| Nucleus maximum eccentricity | 0.80 | 0.85 | V2 configuration of 14 January; 0.80 removes 31 of the original cells |
| Nucleus proximity exclusion | on | off | not present in the code of January 2026 |
| Sato ridge polarity | scikit-image default (dark ridges) | bright ridges (`sato_black_ridges: false`) | the dark-ridge setting traces the gaps between fibres; the preceding tool (Skeletactin, October 2025) used bright ridges |
| Sato σ | 0.5 to 2.0, 5 steps | 0.7, 0.96, 1.22, 1.48, 1.74, 2.0 | values of January 2026 (V2 configuration, code default, confirmed by the authors) |
| Minimum object size / closing | 60 px / 0 | 3 px / 1 | values printed in the notebook of 17 January and in the reference block of 17 February |
| Additional Vieira smoothing | applied (σ 0.5 despite being disabled) | off | disabled in all settings of January 2026 |

Unchanged: nuclei brightness/contrast 1.4 / 1.3; Cellpose `cpsam` with automatic diameter, flow
threshold 0.4, cell probability threshold 0.0; Gaussian σ 0.7 before the Sato filter; Li threshold with
hysteresis (low fraction 0.4); zones at 0.33 and 0.66 of the nucleus-to-membrane distance; cortex
defined as the outer 5 px.

Tested and rejected: the Cellpose settings of the 17 January notebook (flow 0.25, probability −1) and its
brightness/contrast 1.9 / 1.75, which retained fewer of the original cells and gave intensities 1.3 to
1.5 times too high; the custom 3D Sato filter of the notebooks, which showed no cell-level agreement. A
55th-percentile threshold is kept as a sensitivity analysis (weighted densities ×1.12 instead of ×1.21,
branch counts ×0.72).

Cellpose model. The V3 code requests `pretrained_model='cpsam'` (Cellpose-SAM, Cellpose 4). The exact
Cellpose version of the original runs was not recorded. With Cellpose 3.1, which was installed later,
the name `cpsam` is not recognized and Cellpose falls back to `cyto3` without an error. The model was
therefore determined on the 18 representative Exp10 fields, with all other settings fixed:

| Model | Cells retained | Median cell area (px) | Original cells recovered |
|---|---|---|---|
| Original analysis | 252 | 19,704 | not applicable |
| Cellpose 3, as coded (falls back to `cyto3`) | 303 | 6,153 | 0 |
| `cyto2` cells, `nuclei` nuclei | 97 | 9,330 | 6 |
| `cyto2` cells, `cyto2` nuclei | 118 | 9,128 | 5 |
| Cellpose 4 `cpsam` | 256 | 20,170 | 211 (84 %), cell area r = 0.999 |

Only `cpsam` reproduces the original segmentation, which identifies it as the model used.

`run_v3.py` applies the two code-level changes (intensities on the adjusted image, bright-ridge Sato
filter) at run time through `replication/v3_paper.py`; the files in `old_pipeline/` are not modified.
Each run writes its settings to `settings_used.json`.

## Agreement with the original analysis

Re-runs with `config_paper.json` through `run_v3.py`, compared with the original tables by
`compare_to_reference.py` (2026-10-09). Recovered cells are original cells matched by area and
perimeter in the re-run. Ratios are medians of re-run/original over recovered cells; r is the
cell-level Pearson correlation.

| | Exp10, 18 fields (calibration) | Exp1, 10 fields (held out) | Exp8, 10 fields (held out) |
|---|---|---|---|
| Cells retained, original / re-run | 252 / 256 | 119 / 120 | 218 / 209 |
| Original cells recovered | 211 (84 %) | 110 (92 %) | 167 (77 %) |
| cell_area | 1.00 (r 0.999) | 1.00 (r 0.999) | 1.01 (r 0.999) |
| cortical_intensity | 1.07 (r 0.96) | 1.11 (r 0.99) | 1.06 (r 0.88) |
| inner_intensity | 1.03 (r 0.96) | 1.08 (r 0.99) | 1.08 (r 0.92) |
| cortical_weighted_density | 1.20 (r 0.84) | 0.98 (r 0.89) | 1.15 (r 0.82) |
| skeleton_density | 1.11 | 1.11 | 1.16 |
| fidelity | 1.30 | 1.30 | 1.39 |

With the February settings, the same Exp10 fields gave 183 cells, 56 % recovery, `cortical_intensity`
×0.88, `cortical_weighted_density` ×0.70 and branch counts ×1.5.

Between 77 and 92 % of the original cells are recovered, and the total number of retained cells is
similar. For about 8 % of the original cells no matching Cellpose mask is produced; the exact Cellpose 4
release of the original runs was not recorded. Cortical intensity and cortical weighted density agree
with the original values within 2 to 20 % depending on the plate and correlate at the cell level.
Fidelity is higher in the re-run (×1.3 to 1.4) and does not correlate at the cell level; all re-run
cells pass the fidelity threshold of 0.25, compared with about 97 % originally. The fine structure of
the skeleton is therefore not fully reproduced, although branch counts agree (×1.04, r 0.94, Exp10).
With two fields per held-out well, most well medians are within 10 % of the original and a few differ
by 20 to 40 %. Whole-plate re-runs (section 4) with 10 fields per well are required for a well-level
comparison.

The reported result is reproduced if the re-run plates show the same direction and significance of
BasR − KpSen in the statistical analysis of section 1 (per-plate tests and pooled mixed model) and well
medians in the range of figure panels F and G.

## Troubleshooting

| Message or symptom | Cause and solution |
|---|---|
| `No *_Z_stack_merged.tif(f) files found` | Point `prepare_input.py` to the plate folder or its `Z_Stack_Merged/` subfolder. Raw Harmony `Images/` contain single planes and are not read directly. |
| `WARNING: … listed field(s) are not on this plate` | The field list belongs to another plate. |
| Cellpose runs on the CPU although a GPU is present | The installed PyTorch build has no CUDA support. Reinstall `torch` from pytorch.org for the installed CUDA driver. |
| GPU out of memory | One field requires about 3.4 GB of GPU memory; close other GPU processes. |
| Intensities 10 to 25 times lower than the original | The run used `--measure-on raw`, or `old_pipeline/main_pipeline_v3.py` was called directly. Use `replication/run_v3.py`. |
| All cells assigned to group "Unassigned" | File names do not start with `<Row><Column>_` (for example `C9_`), or the plate layout differs; edit `well_to_group`. |

## Licence

MIT, see `LICENSE`.
