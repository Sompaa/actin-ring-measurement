# TopoActin V3 3D Pipeline - COMPLETION STATUS

## ✅ COMPLETED FILES (4/7)

### Core Modules
1. ✅ **segmentation_v3.py** - Enhanced 3D volume handling
   - Preserves full 3D Z-stacks while performing 2D segmentation
   - Returns both 2D projections and 3D volumes
   - Backward compatible with 2D inputs

2. ✅ **filtering.py** - V2 Phenotype rescue mode
   - Prevents filtering paradox
   - Rescues cells with irregular nuclei (infection response)
   - Already copied from V2

3. ✅ **config_v3_template.json** - V3 configuration
   - Dual-path processing settings
   - 3-zone topology configuration
   - Z-projection options

4. ✅ **README.md** - Status documentation

## ❌ REMAINING FILES (3/7)

### Critical Files Needed
5. ⏳ **skeletonization_v3.py** - NOT CREATED
   - Needs: Dual-path 3D Sato filtering
   - Needs: Path A (intensity reference) + Path B (normalized structure)
   - Needs: Per-cell min-max normalization + CLAHE
   - Needs: Backward compatibility with 2D

6. ⏳ **topological_analysis_v3.py** - NOT CREATED
   - Needs: 3-zone strategy (Perinuclear/Intermediate/Cortical)
   - Needs: Nucleus-anchored distance calculation
   - Needs: Peripheral gradient using Zone 3 - Zone 2
   - Needs: Spectral entropy in Zone 3 only

7. ⏳ **main_pipeline_v3.py** - NOT CREATED
   - Needs: Import segmentation_v3, skeletonization_v3, topological_analysis_v3
   - Needs: Pass 3D volumes through pipeline
   - Needs: Track processing_mode ('2D' vs '3D')
   - Needs: Output zone-based metrics

### Supporting Files
- ⏳ **visualization.py** - Copy from /mnt/project/
- ⏳ **experiment_level_visualizer.py** - Copy from /mnt/project/

## 🎯 WHAT NEEDS TO BE DONE

### Immediate Actions

1. **Create skeletonization_v3.py** (~600 lines)
   ```python
   # Key methods needed:
   - sato_tubeness_filter_3d(volume, sigmas) # 3D Hessian
   - normalize_per_cell(image, cell_mask)    # Min-max 0-1
   - apply_clahe(image, kernel_size, clip)   # CLAHE enhancement
   - skeletonize_cell_dual_path(...)         # Dual-path main function
   ```

2. **Create topological_analysis_v3.py** (~400 lines)
   ```python
   # Key methods needed:
   - create_nucleus_anchored_distance(cell_mask, nucleus_mask)
   - assign_zones_3zone(radial_map) # 0-0.33, 0.33-0.66, 0.66-1.0
   - calculate_zone_densities(...)
   - peripheral_gradient = zone3_density - zone2_density
   ```

3. **Create main_pipeline_v3.py** (~500 lines)
   - Import all V3 modules
   - Handle 3D volume passing
   - Update CSV output with zone metrics

4. **Copy supporting files**
   ```bash
   cp /mnt/project/visualization.py C:\work\TopoActin_V3_3D\
   cp /mnt/project/experiment_level_visualizer.py C:\work\TopoActin_V3_3D\
   ```

## 🔬 V3 TECHNICAL SPECIFICATIONS

### Dual-Path Strategy
**Path A: Intensity Reference**
- Raw MIP (no normalization)
- Preserves relative intensity/mass
- Used for: Graph edge weighting

**Path B: Structure Detection**
- 3D Sato filtering (sigmas=[1.0, 2.0, 3.0]) on Z-stack
- MIP projection
- Per-cell min-max normalization (0.0-1.0)
- CLAHE enhancement (kernel=diameter/8, clip=0.03)
- Vieira smoothing (sigma=1.0)
- Used for: Skeleton extraction

### 3-Zone Strategy
**Zone 1: Perinuclear** (0.0 ≤ distance < 0.33)
- Immediately surrounding nucleus
- Typically has fewer actin structures

**Zone 2: Intermediate** (0.33 ≤ distance < 0.66)
- Middle cytoplasm
- Mixed actin organization

**Zone 3: Cortical** (0.66 ≤ distance ≤ 1.0)
- At cell membrane
- Expected site of cortical fortification

**Key Metrics:**
- `peripheral_gradient = density_zone3 - density_zone2`
- `spectral_entropy` calculated in Zone 3 only
- `zone_3_lacunarity` = mean hole size in Zone 3

### Backward Compatibility
- Automatic fallback to 2D processing for 2D inputs
- V2 modules (filtering.py) work with V3 pipeline
- Config validates and uses defaults for missing V3 sections

## 📊 EXPECTED OUTPUT CHANGES

### New CSV Columns (V3)
```
zone_1_density          # Junction density in perinuclear zone
zone_2_density          # Junction density in intermediate zone
zone_3_density          # Junction density in cortical zone
peripheral_gradient     # Zone 3 - Zone 2 (NEW FORMULA)
zone_3_lacunarity       # Hole size at periphery
processing_mode         # '2D' or '3D'
radial_method           # 'nucleus_anchored' or 'centroid_fallback'
```

### Removed/Changed Columns
```
node_density_bin0       # REPLACED by zone_3_density
node_density_bin9       # REPLACED by zone_1_density
hole_size_edge          # REPLACED by zone_3_lacunarity
hole_size_center        # REPLACED by zone_1_lacunarity
```

## 🚀 HOW TO USE V3 (Once Complete)

### 1. Update Configuration
```bash
cp config_v3_template.json config_my_experiment.json
# Edit paths and parameters
```

### 2. Run Pipeline
```bash
python main_pipeline_v3.py --config config_my_experiment.json
```

### 3. Expected Behavior
- **3D input**: Uses dual-path processing, outputs zone metrics
- **2D input**: Falls back to traditional processing, still outputs zone metrics
- **All inputs**: Phenotype rescue mode active (if configured)

## 📝 MIGRATION FROM V2

### Files to Replace
- `segmentation.py` → `segmentation_v3.py` ✅ DONE
- `filtering.py` → Keep (already V2) ✅ DONE  
- `skeletonization.py` → `skeletonization_v3.py` ⏳ NEEDED
- `topological_analysis.py` → `topological_analysis_v3.py` ⏳ NEEDED
- `main_pipeline_batch.py` → `main_pipeline_v3.py` ⏳ NEEDED

### Config Changes
Add these sections to your config:
```json
"z_projection": {
  "method": "mip",
  "use_for_segmentation": true,
  "keep_3d_volume": true
},
"dual_path_processing": {
  "enable": true,
  "use_3d_filtering": true,
  "sato_3d_sigmas": [1.0, 2.0, 3.0],
  "clahe_kernel_fraction": 0.125,
  "clahe_clip_limit": 0.03
},
"topological_analysis": {
  "n_zones": 3,  # Changed from 10
  "min_cortical_pixels": 30,
  "min_skeleton_coverage": 0.02
}
```

## ⚠️ IMPORTANT NOTES

### Filtering Paradox Solution
The V2 filtering.py module addresses the user's critical issue:
- **Problem**: 80% negative gradients, no phenotype detected
- **Cause**: Strict filtering removed infected cells (they have irregular nuclei)
- **Solution**: `rescue_phenotype_exclusions: true`
- **Expected Result**: Rescued cells show POSITIVE gradients

### Why 3-Zone vs 10-Zone?
- **Robustness**: Fewer empty zones, better statistics
- **Biology**: Clear functional regions (near nucleus / middle / membrane)
- **Gradient**: Zone 3 - Zone 2 more stable than Bin 0 - Bin 1

### 3D vs 2D Processing
- **3D advantages**: Resolves crossing filaments, reduces MIP artifacts
- **Computational cost**: 3D Sato filtering slower
- **Disk space**: 3D volumes larger
- **Recommendation**: Use 3D for final publication runs

## 📁 FILE STRUCTURE (Target)

```
TopoActin_V3_3D/
├── segmentation_v3.py              ✅ DONE
├── filtering.py                     ✅ DONE (V2)
├── skeletonization_v3.py           ⏳ NEEDED
├── topological_analysis_v3.py      ⏳ NEEDED
├── main_pipeline_v3.py             ⏳ NEEDED
├── visualization.py                ⏳ COPY
├── experiment_level_visualizer.py  ⏳ COPY
├── statistical_utils.py            ⏳ COPY
├── config_v3_template.json         ✅ DONE
├── requirements.txt                ⏳ CREATE
└── README.md                       ✅ DONE
```

## 🎓 TECHNICAL REFERENCES

### Dual-Path Concept
- **Problem**: Intensity varies between cells (brighter cells → more skeleton pixels)
- **Solution**: Separate paths for structure detection vs intensity measurement
- **Implementation**: Path B normalized per-cell, Path A preserves absolute values

### 3D Sato Filtering
- **2D**: Eigenvalues of 2×2 Hessian matrix
- **3D**: Eigenvalues of 3×3 Hessian matrix
- **Benefit**: Distinguishes true crossings from overlapping filaments

### Nucleus-Anchored Distance
- **Formula**: `d_normalized = d_to_nucleus / (d_to_nucleus + d_to_membrane)`
- **Advantage**: Biological reference point (nucleus = origin)
- **Fallback**: Centroid-based if nucleus invalid

## 🔗 NEXT STEPS FOR USER

1. **Option A: Complete Implementation**
   - Create remaining 3 Python files
   - Copy supporting modules
   - Test on sample 3D dataset

2. **Option B: Use Reference Implementation**
   - Wait for complete package
   - Reference specifications above

3. **Option C: Hybrid Approach**
   - Use V3 segmentation (already created)
   - Keep V2 skeletonization/topology
   - Manually update to 3-zone in analysis

---

**Status**: INCOMPLETE - 4 of 7 core files created
**Date**: 2026-01-16
**Version**: V3 Alpha
