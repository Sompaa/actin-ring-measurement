# Implementation Progress: Bacterial Burden Integration

## Plan Reference

See full plan at: `C:\Users\arons\.claude\plans\wobbly-leaping-rose.md`

## Status Overview

| Section | Description | Status | Files Changed |
| --- | --- | --- | --- |
| 7 | IMPLEMENTATION_PROGRESS.md | COMPLETE | this file |
| 1 | ARCHITECTURE.md | COMPLETE | ARCHITECTURE.md (new) |
| 2 | Config schema extension | COMPLETE | config_v3_template.json |
| 3 | Extract detection functions | COMPLETE | bacterial_burden/detection.py (new), segmentation_bacterial.py |
| 4 | Fix Cellpose loading bug | COMPLETE | segmentation_v3.py, main_pipeline_v3.py |
| 5 | Extend segmentation for N-channel | COMPLETE | segmentation_v3.py (full rewrite) |
| 6 | Main pipeline bacterial integration | COMPLETE | main_pipeline_v3.py |
| EXTRA | Interactive notebook preview | COMPLETE | TopoActin_V3_Interactive_DOWNSAMPLING.ipynb |

## Dependency Order

```text
Section 7 (this file)         -- DONE
Section 1 (ARCHITECTURE.md)   -- DONE
Section 2 (config schema)     -- DONE
Section 3 (detection.py)      -- DONE
Section 4 (Cellpose fix)      -- DONE
Section 5 (segmentation ext)  -- DONE (depends on 4)
Section 6 (main pipeline)     -- DONE (depends on 3 + 5)
```

## Change Log

### Section 7 - IMPLEMENTATION_PROGRESS.md

- Created this tracking file

### Section 1 - ARCHITECTURE.md

- Created comprehensive architecture reference document
- Covers: file inventory, data flow, channel handling, config schema, API summaries
- Documents 16 known vulnerabilities/bad practices (5 fixed, 11 still present)
- Includes experimental design context (plate layout, strains, key metrics)

### Section 2 - Config Schema Extension

- Added `bacterial_burden` block to `config_v3_template.json`
- Fields: enabled, channels (actin/nuclei/bacteria/brightfield indices), cell_segmentation_channel, brightness_contrast (bacteria/brightfield), bacteria_detection (sato params), per_cell_normalization, qc
- Default: disabled (backward compatible)
- Per-experiment overrides supported in batch mode

### Section 3 - Extract Detection Functions

- Created `bacterial_burden/detection.py` with 3 standalone functions:
  - `normalize_cell_region()` - per-cell intensity normalization
  - `detect_bacteria_in_cell()` - Sato ridge filter on bacteria ROI
  - `quantify_bacteria_per_cell()` - orchestrator loop over cells
- Updated `segmentation_bacterial.py`:
  - `_normalize_cell_region()` now delegates to `detection.normalize_cell_region()`
  - `_detect_bacteria_sato()` now delegates to `detection.detect_bacteria_in_cell()`
  - `_quantify_bacteria_per_cell()` unchanged (uses self methods which now delegate)

### Section 4 - Fix Cellpose Loading Bug

- `segmentation_v3.py`: Moved `CellposeModel` creation from `segment_image()` to `__init__()`
  - `self.cell_model` and `self.nucleus_model` loaded once
  - `self.cell_cfg` and `self.nucleus_cfg` stored for eval calls
- `main_pipeline_v3.py`: Moved `CellSegmenter(config)` from inside per-image loop to before loop

### Section 5 - Extend segmentation_v3.py

- Full rewrite of `segmentation_v3.py` (246 -> ~280 lines)
- `__init__()`: Reads channel indices from `bacterial_burden.channels` config
  - `ch_actin`, `ch_nuclei`, `ch_bacteria`, `ch_brightfield`, `cell_seg_channel`
  - Backward compatible: defaults to actin=0, nuclei=1 when bacterial_burden disabled
- New `_validate_image()`: Checks empty, ndim, dtype, all-zero
- New `_validate_channel_count()`: Ensures image has enough channels for config
- `load_and_preprocess()`:
  - Now returns dict instead of 7-tuple
  - Config-driven channel extraction (no more hardcoded `[:, 0]` / `[:, 1]`)
  - Extracts bacteria_2d via MIP when enabled
  - Extracts brightfield_2d when configured
  - Applies B/C to bacteria and brightfield channels
- `segment_image()`:
  - Uses brightfield for Cellpose when `cell_segmentation_channel == "brightfield"`
  - Falls back to actin if brightfield not available
  - Returns bacteria_2d, bacteria_2d_raw, brightfield_2d in result dict

### Section 6 - Main Pipeline Bacterial Integration

- Added conditional import of `bacterial_burden.detection.quantify_bacteria_per_cell`
- `process_experiment()`:
  - Builds `effective_config` with per-experiment bacterial_burden overrides
  - Determines `bacterial_enabled` flag
  - After filtering step: runs `quantify_bacteria_per_cell()` on bacteria channel
  - Inside per-cell loop: merges `intracellular_bact_area`, `area_ratio`, `has_bacteria` into row_data
  - When bacterial_burden disabled: no extra CSV columns

---

## Verification Checklist

- [ ] Backward compat: Run with existing config (no bacterial_burden key) on 2-channel data
- [ ] 3-channel bacterial: Run with bacterial_burden.enabled=true on 3-channel MIP
- [ ] Brightfield mode: Run with cell_segmentation_channel="brightfield" on 4-channel data
- [ ] Batch mode: Mixed experiments with different bacterial_burden settings
- [ ] Standalone: main_bacterial_burden.py still works after detection.py refactor
- [ ] Performance: Cellpose loads once per experiment

---

### EXTRA - Interactive Notebook Preview

- **cell-7 (Image Loading)**: Updated to detect and extract bacteria (ch2) and brightfield (ch3) from 3/4-channel images, displays all available channels in dynamic panel layout
- **New markdown cell after cell-14**: "Step 5b: Bacterial Burden Preview (Optional)" header
- **New code cell after markdown**: Full bacterial burden preview with interactive widgets for Sato parameters (sato_sigmas, threshold_method, min_area, min_signal_range, padding), runs `quantify_bacteria_per_cell()`, shows RGB overlay (Red=bacteria, Green=actin, Blue=cell boundary) for top 6 infected cells
- **cell-28 (Config Export)**: Added `bacterial_burden` section to exported config JSON, auto-detects bacteria/brightfield availability, captures widget parameter values from Step 5b

---

## Post-Implementation Changes

### Channel Index Swap (2026-02-16)

User requested swapping brightfield and bacteria channel indices throughout the codebase:

**Old order (4-channel):** actin(0), nuclei(1), bacteria(2), brightfield(3)

**New order (4-channel):** actin(0), nuclei(1), **brightfield(2)**, **bacteria(3)**

**Files updated:**

- `config_v3_template.json`: channels block (line 165-168)
- `TopoActin_V3_Interactive_DOWNSAMPLING.ipynb`:
  - cell-7: 4-channel extraction logic (brightfield=ch2, bacteria=ch3)
  - cell-28: config export channel assignments
- `segmentation_v3.py`: No changes needed (reads from config)

**Backward compatibility:** 3-channel images (actin + nuclei + bacteria) unchanged - bacteria remains channel 2.

---

*If context resets, read this file first, then ARCHITECTURE.md, then continue from the next pending item.*
