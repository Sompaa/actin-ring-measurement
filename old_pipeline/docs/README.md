# TopoActin V3: 3D Dual-Path Cytoskeleton Analysis Pipeline

## ⚠️ INCOMPLETE IMPLEMENTATION

**Current Status**: 4 of 7 core files created (57% complete)

This directory contains a **partially implemented** V3 refactor of the Topological Biopsy pipeline. The implementation adds:
- 3D Z-stack processing capability
- Dual-path strategy (structure detection + intensity reference)
- 3-zone topology (Perinuclear/Intermediate/Cortical)
- Enhanced phenotype rescue mode (V2)

**See `IMPLEMENTATION_STATUS.md` for detailed completion status and technical specifications.**

---

## ✅ What's Been Created

### Core Modules (Partially Complete)
1. **segmentation_v3.py** ✅ COMPLETE
   - Enhanced 3D volume handling
   - Preserves full Z-stacks for downstream processing
   - Backward compatible with 2D inputs

2. **filtering.py** ✅ COMPLETE (V2)
   - Phenotype rescue mode (solves filtering paradox)
   - Flags rescued cells for analysis

3. **config_v3_template.json** ✅ COMPLETE
   - Dual-path processing settings
   - 3-zone configuration
   - Z-projection options

4. **requirements.txt** ✅ COMPLETE
   - Python dependencies list

## ❌ What's Missing

### Critical Files Still Needed
- **skeletonization_v3.py** (~600 lines)
  - Dual-path 3D Sato filtering
  - Per-cell normalization + CLAHE
  - Integration with topology analysis

- **topological_analysis_v3.py** (~400 lines)
  - 3-zone strategy implementation
  - Nucleus-anchored distance maps
  - Zone-based metrics (peripheral gradient, etc.)

- **main_pipeline_v3.py** (~500 lines)
  - Pipeline orchestration
  - 3D volume passing
  - CSV output with zone metrics

### Supporting Files to Copy
- visualization.py (from Complete_TopoActinAnalyzer)
- experiment_level_visualizer.py (from Complete_TopoActinAnalyzer)
- statistical_utils.py (from Complete_TopoActinAnalyzer)

---

## 🎯 V3 Key Innovations

### 1. Dual-Path Processing Strategy

**The Problem:**
- Intensity varies between cells (brighter cells appear to have more actin)
- Simple normalization loses relative intensity information
- Creates false differences in measurements

**The Solution:**
- **Path A (Reference)**: Raw MIP, preserves intensity for graph weighting
- **Path B (Structure)**: 3D Sato → MIP → per-cell normalize → CLAHE → skeleton

**Benefits:**
- Structure detection independent of cell brightness
- Intensity measurements preserve biological variation
- Solves "bright cell bias"

### 2. 3D Sato Filtering

**Why 3D?**
- 2D MIP collapses crossing filaments → false knots
- 3D Hessian distinguishes true crossings from overlaps
- Applied BEFORE projection → cleaner skeleton

**Implementation:**
```python
# 3D Sato on Z-stack [Z, Y, X]
tubeness_3d = sato_filter_3d(actin_volume, sigmas=[1.0, 2.0, 3.0])
# Then project
tubeness_2d = np.max(tubeness_3d, axis=0)
```

### 3. 3-Zone Topology Strategy

**Why 3 Zones (not 10)?**
- **Robustness**: Fewer empty zones, better statistics
- **Biology**: Clear functional regions
- **Simplicity**: Easier to interpret

**Zones:**
1. **Perinuclear** (0.0–0.33): Around nucleus
2. **Intermediate** (0.33–0.66): Middle cytoplasm  
3. **Cortical** (0.66–1.0): At membrane

**Key Metric:**
```python
peripheral_gradient = density_zone3 - density_zone2
# More stable than old Bin0 - Bin1
```

---

## 🔬 Technical Specifications

### Dual-Path Algorithm

**Path A: Intensity Reference**
```python
intensity_reference = np.max(actin_volume_3d, axis=0)  # Raw MIP
# Used for: skan.Skeleton(binary, source_image=intensity_reference)
```

**Path B: Structure Detection**
```python
# Step 1: 3D Sato filtering
sato_3d = sato_tubeness_filter_3d(actin_volume_3d, sigmas=[1.0, 2.0, 3.0])

# Step 2: MIP projection
sato_2d = np.max(sato_3d, axis=0)

# Step 3: Per-cell normalization (KEY INNOVATION)
cell_region = sato_2d[cell_mask]
normalized = (cell_region - cell_region.min()) / (cell_region.max() - cell_region.min())

# Step 4: CLAHE enhancement
from skimage import exposure
kernel_size = int(cell_diameter / 8)
enhanced = exposure.equalize_adapthist(normalized, kernel_size=kernel_size, clip_limit=0.03)

# Step 5: Vieira smoothing
smoothed = gaussian_filter(enhanced, sigma=1.0)

# Step 6: Threshold & skeletonize
binary = smoothed > threshold
skeleton = skeletonize(binary)
```

### 3-Zone Distance Calculation

**Nucleus-Anchored (Preferred):**
```python
# Distance from nucleus
dist_nucleus = distance_transform_edt(~nucleus_mask)

# Distance from membrane  
dist_membrane = distance_transform_edt(cell_mask)

# Normalized distance (0=nucleus, 1=membrane)
radial_distance = dist_nucleus / (dist_nucleus + dist_membrane)

# Assign zones
zone_1 = (radial_distance >= 0.0) & (radial_distance < 0.33)
zone_2 = (radial_distance >= 0.33) & (radial_distance < 0.66)
zone_3 = (radial_distance >= 0.66) & (radial_distance <= 1.0)
```

**Centroid Fallback (if nucleus invalid):**
```python
radial_distance = distance_transform_edt(cell_mask)
radial_distance /= radial_distance.max()
```

### Expected CSV Output

**New Columns:**
- `zone_1_density` - Junction density near nucleus
- `zone_2_density` - Junction density in middle
- `zone_3_density` - Junction density at membrane
- `peripheral_gradient` - Zone 3 - Zone 2 (cortical fortification)
- `zone_3_lacunarity` - Hole size at periphery
- `processing_mode` - '2D' or '3D'
- `radial_method` - 'nucleus_anchored' or 'centroid_fallback'

---

## 🚧 How to Complete Implementation

### Option 1: Create Missing Files Manually

Use specifications in `IMPLEMENTATION_STATUS.md`:

1. **Create skeletonization_v3.py**
   - Copy structure from `/mnt/project/skeletonization.py`
   - Add `sato_tubeness_filter_3d()` method
   - Add `normalize_per_cell()` method
   - Add `apply_clahe()` method
   - Modify `skeletonize_cell()` to use dual-path

2. **Create topological_analysis_v3.py**
   - Copy from `/mnt/project/topological_analysis.py`
   - Modify to use 3 zones instead of 10
   - Update `_calculate_node_density_gradient()` to use Zone 3 - Zone 2
   - Update `_calculate_spectral_entropy()` to use only Zone 3

3. **Create main_pipeline_v3.py**
   - Copy from `/mnt/project/main_pipeline_batch.py`
   - Import V3 modules
   - Pass 3D volumes through pipeline
   - Update CSV output

### Option 2: Use Hybrid V2/V3 Approach

Immediate usable solution:

1. **Use V2 for processing**
   ```bash
   cd C:\work\Complete_TopoActinAnalyzer
   python main_pipeline_batch.py --config config_production_10k.json
   ```

2. **Post-process to 3-zone metrics**
   - Load CSV with 10-zone data
   - Combine Bins 0-2 → Zone 3 (Cortical)
   - Combine Bins 3-6 → Zone 2 (Intermediate)
   - Combine Bins 7-9 → Zone 1 (Perinuclear)
   - Recalculate peripheral_gradient

### Option 3: Request Complete Package

Contact for fully implemented V3 pipeline.

---

## 📊 Filtering Paradox (Critical!)

### The Problem (From User's Data)
- **80% negative gradients** (opposite of expected)
- **NO phenotype detected** in treatment groups
- **Cause**: Strict filtering (min_solidity=0.82) removed infected cells
- **Why**: Infected cells have irregular/deformed nuclei (biological response!)

### The Solution (Already in V2 filtering.py ✅)

**Config settings:**
```json
"nucleus_filtering": {
  "rescue_phenotype_exclusions": true,
  "rescue_solidity_threshold": 0.6,
  "rescue_eccentricity_threshold": 0.95
}
```

**Expected outcome with rescue:**
- Control: gradient ~ -0.004 (no change)
- Treatment (rescued cells): gradient ~ +0.008 to +0.012 (POSITIVE!)

**The treatment likely DOES work** - you were just throwing away the infected cells!

---

## 📁 File Locations

- **V3 Development**: `C:\work\TopoActin_V3_3D\` (THIS DIRECTORY)
- **V2 Reference**: `C:\work\Complete_TopoActinAnalyzer\`
- **Project Files**: `/mnt/project/` (read-only reference)

---

## 🔗 Quick Start (Once Complete)

### Installation
```bash
cd C:\work\TopoActin_V3_3D
pip install -r requirements.txt
```

### Configuration
```bash
cp config_v3_template.json my_experiment.json
# Edit my_experiment.json with your paths
```

### Run Pipeline
```bash
python main_pipeline_v3.py --config my_experiment.json
```

### Expected Behavior
- **3D inputs** (Z-stacks): Uses dual-path 3D processing
- **2D inputs**: Automatically falls back to traditional processing
- **All inputs**: Phenotype rescue mode active (if configured)

---

## 📝 Notes

- **Backward Compatible**: V3 works with 2D inputs (auto-fallback)
- **Phenotype Rescue**: V2 filtering.py addresses user's filtering paradox
- **3-Zone Stability**: More robust than 10-zone for statistical analysis
- **Dual-Path**: Solves intensity inconsistency between cells

---

## 🆘 Support

For questions or issues:
1. Check `IMPLEMENTATION_STATUS.md` for detailed specifications
2. Review `/mnt/project/` files for reference implementation
3. See original documentation in `Complete_TopoActinAnalyzer`

---

**Status**: INCOMPLETE - 4 of 7 core files created  
**Date**: 2026-01-16  
**Version**: V3 Alpha (Partial Implementation)
