# Example fields

Four fields of Exp10 (2026-01-13, plate 4; A549, 1 h infection), the fields of the paper figure slides.
Each is a merged Z-stack as V3 reads it: 5 Z-planes × 2 channels (0 = actin, CellMask; 1 = nuclei,
Hoechst), 1080 × 1080 px, uint16.

| File | Group | In the January analysis |
|---|---|---|
| `C9_F13_Z_stack_merged.tiff` | BasR (basR G53V) | yes, 18 cells; paper figure cell = cell 21 in the walkthrough |
| `D10_F02_Z_stack_merged.tiff` | BasR (basR G53V) | no (shown in figures, not in the 10 analysed fields of its well) |
| `B5_F11_Z_stack_merged.tiff` | KpSen (ancestor) | yes, 16 cells |
| `D11_F09_Z_stack_merged.tiff` | Control (uninfected) | yes, 21 cells |

`reference_january_Exp10.csv` holds the January 2026 result rows of these fields (from the original
Exp10 table), to compare a re-run with (`replication/compare_to_reference.py`).

Quick test of the whole procedure on one field:

```bash
jupyter lab notebooks/single_image_walkthrough.ipynb     # uses examples/C9_F13 by default
# or the batch pipeline on all four:
python replication/run_v3.py examples/ out_examples/ --rings
python replication/compare_to_reference.py out_examples/results/cytoskeleton_analysis_v3_*.csv examples/reference_january_Exp10.csv -o compare_examples/
```

For C9_F13 cell 21 the walkthrough gives inner / mid / outer Sato-weighted density 0.0015 / 0.0028 /
0.0048, the values printed on the February 2026 figure slide for that cell.
