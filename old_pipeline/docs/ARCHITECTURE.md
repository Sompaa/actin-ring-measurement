# TopoActin V3 Pipeline - Architecture Reference

## Project Overview

TopoActin V3 is a cytoskeleton analysis pipeline for fluorescence microscopy images of human cell lines. It quantifies actin network topology using segmentation, skeletonization, and zonal analysis. A sub-module handles bacterial burden quantification for host-pathogen interaction studies.

## File Inventory

### Core Pipeline (root directory)

- `main_pipeline_v3.py` (377 lines) - Pipeline orchestrator. Batch/single mode, CSV export, visualization
- `segmentation_v3.py` (246 lines) - CellSegmenter class. TIFF loading, Z-projection, Cellpose segmentation, B/C/gamma adjustment
- `filtering.py` (378 lines) - NucleiFilter + CellFilter classes. Morphology filters, V2 phenotype rescue mode
- `skeletonization_v3.py` (630 lines) - CellSkeletonizer class. Dual-path Sato processing, weighted skeleton metrics
- `topological_analysis_v3.py` (193 lines) - TopologicalAnalyzer class. 3-zone nucleus-anchored topology

### Visualization

- `visualization.py` (199 lines) - ResultsVisualizer. Per-cell comparison PNG
- `experiment_level_visualizer.py` (387 lines) - ExperimentLevelVisualizer. Violin/stats plots for 1000+ cells
- `create_pptx.py` (730 lines) - PowerPoint data flow generator

### Utilities

- `position_downsampling.py` (424 lines) - Well-based position sampling with reproducible seeds

### Bacterial Burden Module (`bacterial_burden/`)

- `main_bacterial_burden.py` (403 lines) - Standalone bacterial pipeline orchestrator
- `segmentation_bacterial.py` (648 lines) - BacterialSegmenter. 3-channel processing, Sato bacteria detection
- `filtering_bacterial.py` (204 lines) - HostCellFilter. 4-stage hierarchical filtering
- `quantification_bacterial.py` (84 lines) - BacterialBurdenQuantifier. Per-cell metrics assembly
- `plotting.py` (364 lines) - Violin/box plots for bacterial burden results
- `detection.py` (NEW - being created) - Extracted standalone detection functions

### Config

- `config_v3_template.json` - V3 pipeline configuration template
- `config_v3.json` - Active experiment config
- `bacterial_burden/config_bacterial_template.json` - Bacterial module template
- `bacterial_burden/config_bacterial.json` - Active bacterial config

### Interactive Notebooks

- `TopoActin_V3_Interactive_DOWNSAMPLING.ipynb` - Interactive parameter tuning with downsampling
- `INTERACTIVE_PREVIEW_COMPLETE.ipynb` - Full interactive preview
- `bacterial_burden/BacterialBurden_Interactive_Preview.ipynb` - Bacterial parameter tuning

## Pipeline Data Flow

```
Image File (.tif/.tiff)
    |
    v
[segmentation_v3.py] CellSegmenter.segment_image()
    |- Load TIFF (2D/3D/4D auto-detect)
    |- Z-projection (MIP/mean/single_middle)
    |- Brightness/contrast/gamma adjustment
    |- Cellpose cpsam: cell masks + nucleus masks
    |
    v
[filtering.py] NucleiFilter.filter_nuclei()
    |- Area bounds, solidity, eccentricity, roundness
    |- V2 Phenotype Rescue: recovers deformed nuclei (biological, not artifact)
    |
    v
[filtering.py] CellFilter.filter_cells()
    |- Area bounds, solidity, eccentricity
    |- Border exclusion, nucleus requirement
    |
    v
[skeletonization_v3.py] CellSkeletonizer.skeletonize_cell() (per cell)
    |- Path A (Intensity): Raw MIP -> intensity reference
    |- Path B (Structure): Sato tubeness -> normalize -> CLAHE -> threshold -> skeleton
    |- Weighted skeleton density in 3 radial zones
    |
    v
[topological_analysis_v3.py] TopologicalAnalyzer.analyze_topology() (per cell)
    |- Nucleus-anchored distance map (0=nucleus, 1=membrane)
    |- 3 zones: perinuclear [0-0.33], intermediate [0.33-0.66], cortical [0.66-1.0]
    |- Peripheral gradient, lacunarity, spectral entropy
    |
    v
[main_pipeline_v3.py] -> CSV with ~45 columns per cell
    |- experiment, image_name, cell_id, group
    |- Morphology: cell_area, nucleus_area, solidity, eccentricity
    |- Skeleton: density, branches, junctions, fidelity, tortuosity
    |- Weighted: total_network_intensity, mean_bundle_intensity, weighted_total_length
    |- Zonal: zone_1/2/3_density, peripheral_gradient, cortical ratios
    |- Topology: spectral_entropy, zone_3_lacunarity
```

## Channel Handling

| Format | Shape | Channels | Handling |
|--------|-------|----------|----------|
| 2D single | [Y, X] | Actin only | No nuclei, no bacteria |
| 2D multi-channel | [C, Y, X] | C=2: [actin, nuclei] | Direct split |
| Z-stack single | [Z, Y, X] | Actin only | MIP projection |
| Z-stack multi | [Z, C, Y, X] | C=2: [actin, nuclei] | MIP per channel |
| Z-stack 3ch | [Z, C, Y, X] | C=3: [actin, nuclei, bacteria] | MIP + bacterial burden |
| Z-stack 4ch | [Z, C, Y, X] | C=4: [actin, nuclei, bacteria, brightfield] | MIP + brightfield seg |

Channel indices are configurable via `bacterial_burden.channels` in config.

## Config Schema Reference

### Main Pipeline (config_v3_template.json)

Key sections:

- `batch_mode` (bool): Single or multi-experiment
- `input_dir` / `output_dir`: File paths
- `downsampling`: { enabled, num_positions, seed }
- `well_to_group`: Maps well column -> treatment name
- `group_colors`: Treatment -> color for plots
- `brightness_contrast`: { actin: {brightness, contrast, gamma}, nuclei: {...} }
- `segmentation`: { cell: {diameter, flow_threshold, cellprob_threshold, min_size}, nucleus: {...} }
- `nucleus_filtering`: { min/max_area, min_solidity, max_eccentricity, rescue_phenotype_exclusions, rescue thresholds }
- `cell_filtering`: { min/max_area, min_solidity, max_eccentricity, border_buffer, require_valid_nucleus }
- `skeletonization`: { sato_sigmas, threshold_percentile, min_object_size, closing_radius, ... }
- `dual_path_processing`: { enable, use_3d_filtering }
- `topological_analysis`: { n_zones, inner_threshold, mid_threshold }
- `bacterial_burden`: { enabled, channels, cell_segmentation_channel, bacteria_detection, ... }
- `processing`: { save_masks, save_skeletons, use_gpu, batch_size }

### Bacterial Module (config_bacterial_template.json)

Key sections:

- `plate_layout`: Well column -> treatment (same concept as well_to_group)
- `channels`: { host_cytoplasm: 0, nuclei: 1, bacteria: 2 }
- `z_projection`: { method: "best_focus" | "mip" | "single_middle" }
- `qc`: { min_intensity_threshold, saturation_fraction_max }
- `brightness_contrast`: { host, nuclei, bacteria } each with brightness/contrast/gamma
- `segmentation`: { host_cell (cyto3), nucleus, bacteria (bact_fluor_omni) }
- `filtering`: { clear_border, host_cell morphology, require_single_nucleus }
- `bacteria_detection`: { sato_sigmas, threshold_method, min_area, min_signal_range, padding }
- `per_cell_normalization`: { enabled, method, percentile_low, percentile_high }

## Module API Summaries

### CellSegmenter (segmentation_v3.py)

- `__init__(config)` - Loads Cellpose models, reads channel config
- `adjust_brightness_contrast(image, params)` -> uint16
- `apply_z_projection(z_stack, axis)` -> 2D array
- `load_and_preprocess(image_path)` -> dict with 2D projections
- `segment_image(image_path)` -> dict with masks + images

### NucleiFilter / CellFilter (filtering.py)

- `filter_nuclei(nucleus_masks, nuclei_image, image_name)` -> { valid_nuclei, excluded_nuclei, rescued_nuclei }
- `filter_cells(cell_masks, nucleus_masks, valid_nuclei, image_name)` -> { valid_cells, excluded_cells }

### CellSkeletonizer (skeletonization_v3.py)

- `skeletonize_cell(cell_mask, actin_image, ...)` -> { skeleton, features, graph_summary, cortical_features, topological_features, processing_mode }

### TopologicalAnalyzer (topological_analysis_v3.py)

- `get_normalized_distance_map(cell_mask, nucleus_mask)` -> distance map [0,1]
- `analyze_topology(skeleton, distance_map, sato_map)` -> zonal metrics dict

### BacterialSegmenter (segmentation_bacterial.py)

- `load_and_preprocess(image_path)` -> dict with 3-channel 2D + QC
- `segment_host_cells(host_2d, nuclei_2d)` -> labeled masks
- `segment_nuclei(nuclei_2d)` -> labeled masks
- `segment_image(image_path)` -> full result dict with masks + bacteria detection

## Known Vulnerabilities and Bad Practices

### FIXED by Current Implementation Work

1. **Cellpose model loaded per-image** (segmentation_v3.py:187,207) - Models should be loaded once in `__init__`
2. **Hardcoded 2-channel assumption** (segmentation_v3.py:104,121) - `img_stack[:, 0]` and `[:, 1]` instead of config-driven
3. **Fragile 7-tuple return** from `load_and_preprocess` - Changed to dict
4. **No input validation** on image dimensions/dtype/empty files
5. **CellSegmenter created per-image** in main_pipeline_v3.py:213 - Moved outside loop

### Still Present (Technical Debt)

6. **Duplicated B/C adjustment code** - Same algorithm in segmentation_v3.py:37-47 and segmentation_bacterial.py:49-62. Should be a shared utility.
7. **Global variables in interactive notebook** - Heavy reliance on global state for inter-cell communication. Not present in batch pipeline.
8. **Bare `except:` clauses** - main_pipeline_v3.py:306 catches all exceptions silently. Should catch specific exceptions.
9. **Hardcoded paths** - `C:\work\Exp10_Images` in interactive notebook. Should be purely config-driven.
10. **Unicode/emoji in print statements** - Can cause encoding errors on Windows with cp1250 locale.
11. **Mutable default arguments** - `boundaries=[0.33, 0.66]` in notebook helper functions.
12. **No pixel-to-micron calibration** - All area/length measurements are in pixels. No metadata extraction for physical units.
13. **Division by zero risks** - Partially guarded with `+ 1e-10` but not consistently across all normalization functions.
14. **Silent fallback chains** - Import failures (visualization, experiment_level_visualizer) print warnings but execution continues with potentially unexpected behavior.
15. **Memory management** - Matplotlib figures not explicitly closed. `gc.collect()` used sparingly. Large arrays kept in scope.
16. **Cellpose model loaded per-image in bacterial module** - segmentation_bacterial.py:516,548 creates new model instances per call to `segment_host_cells`/`segment_nuclei`.

## Experimental Design Context

### Plate Layout (96-well format)

- Columns 1-12, Rows A-H
- Typical layout: Empty | Control | Strain1 x4 | Strain2 x4 | Control | Empty
- Each well has multiple field positions (F01-F26)
- Well parsing: regex `([A-H])(\d{1,2})_F\d+` from filename

### Bacterial Strains Studied

- **KpSen**: Klebsiella pneumoniae (sensitive strain)
- **BasR**: Klebsiella pneumoniae (resistant strain, BasR mutation)
- Time points: 1h, 3h post-infection
- Replicates: 4 wells per strain

### Key Scientific Metrics

- `peripheral_gradient` = zone_3_density - zone_2_density (cortical fortification)
- `cortical_inner_ratio_weighted` = zone_3 / zone_1 (actin redistribution)
- `area_ratio` = intracellular_bact_area / host_area (bacterial burden per cell)

### V3 Key Innovations

1. **Dual-path processing**: Separates structure detection (Sato, per-cell normalized) from intensity measurement (raw MIP). Eliminates "bright cell bias."
2. **Phenotype rescue mode**: Infected cells have deformed nuclei (biological response). V2 filtering removes them. Rescue mode recovers these strongest phenotypes.
3. **3-zone topology**: Simplified from original 10-zone to 3 zones (perinuclear/intermediate/cortical). More robust statistically, maps to cellular compartments.
4. **Weighted skeleton metrics**: Intensity-weighted skeleton density per zone captures both structural and intensity information.
