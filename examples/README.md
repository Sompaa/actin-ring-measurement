# Example fields

Four fields of Exp10 (2026-01-13, plate 4; A549 cells, 1 h infection), the fields used for the
representative images. Each file is a merged Z-stack in the format read by V3: 5 Z-planes, 2 channels
(0: actin, CellMask; 1: nuclei, Hoechst), 1080 × 1080 px, 16-bit.

| File | Group | Original analysis |
|---|---|---|
| `C9_F13_Z_stack_merged.tiff` | BasR (basR G53V) | analysed, 18 cells; cell 21 in the walkthrough is the cell shown in the figure |
| `D10_F02_Z_stack_merged.tiff` | BasR (basR G53V) | not analysed (illustration only; not among the 10 analysed fields of its well) |
| `B5_F11_Z_stack_merged.tiff` | KpSen (ancestor) | analysed, 16 cells |
| `D11_F09_Z_stack_merged.tiff` | Control (uninfected) | analysed, 21 cells |

`reference_january_Exp10.csv` contains the original per-cell results of these fields, taken from the
Exp10 result table, for comparison with a re-run (`replication/compare_to_reference.py`).

Single-field analysis and figure:

```bash
jupyter lab notebooks/single_image_walkthrough.ipynb     # uses examples/C9_F13 by default
```

Batch pipeline on all four fields, with comparison to the original results:

```bash
python replication/run_v3.py examples/ out_examples/ --rings
python replication/compare_to_reference.py out_examples/results/cytoskeleton_analysis_v3_*.csv examples/reference_january_Exp10.csv -o compare_examples/
```

Expected result: for cell 21 of C9_F13, the walkthrough gives inner, mid and outer Sato-weighted
densities of 0.0015, 0.0028 and 0.0048, the values of the original analysis for this cell.
