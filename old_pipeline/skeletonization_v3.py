"""
Skeletonization Module V3 - Dual-Path 3D Cytoskeleton Processing with Weighted Skeleton Analysis

Key V3 Enhancements:
- Dual-path strategy: Structure detection + Intensity reference
- 3D Sato tubeness filtering (before projection)
- Per-cell normalization (solves brightness inconsistency)
- CLAHE enhancement (adaptive contrast)
- Backward compatible with 2D inputs

NEW V3 Weighted Analysis:
- Uses Sato tubeness (not raw MIP) for weighting skeleton metrics
- Calculates weighted skeleton densities in 3 nucleus-anchored zones
- Computes 3 ratios: Outer/Inner, Outer/Mid, Mid/Inner for comprehensive analysis
- Area-normalized to handle ring size differences

Processing Paths:
  Path A (Reference): Raw MIP -> graph edge weighting (legacy)
  Path B (Structure): 3D Sato -> MIP -> normalize -> CLAHE -> skeleton
  Path C (Weighting): Raw Sato -> skeleton intensity weighting -> zonal analysis
"""

import numpy as np
import tifffile
from pathlib import Path
from scipy import ndimage
from skimage import morphology, filters, exposure
from skimage.morphology import skeletonize, remove_small_objects, binary_closing, disk
from skan import Skeleton, summarize
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Import topological analyzer V3
try:
    from topological_analysis_v3 import TopologicalAnalyzer
    TOPOLOGICAL_AVAILABLE = True
except ImportError:
    TOPOLOGICAL_AVAILABLE = False
    print("Warning: topological_analysis_v3 not found. Topological metrics will be skipped.")


# ============================================================================
# HELPER FUNCTIONS FOR WEIGHTED SKELETON ANALYSIS
# ============================================================================

def calculate_weighted_skeleton_density(skeleton, intensity_map, zone_mask):
    """
    Calculate weighted skeleton density in a zone.
    
    This function computes the total intensity of skeleton pixels within a zone,
    normalized by the zone's total area. This normalization is CRITICAL to account
    for the fact that outer rings have larger area than inner rings.
    
    Args:
        skeleton (ndarray): Binary skeleton mask (Y, X)
        intensity_map (ndarray): Intensity values to weight by (Y, X) - typically Sato tubeness
        zone_mask (ndarray): Boolean mask defining the zone (Y, X)
        
    Returns:
        float: Weighted skeleton density = (sum of intensity on skeleton) / (zone area)
               Returns 0.0 if no skeleton in zone
               Returns np.nan if zone has zero area
               
    Example:
        >>> skeleton = np.array([[1, 0], [0, 1]])  # 2 skeleton pixels
        >>> intensity = np.array([[5, 0], [0, 3]])  # intensities
        >>> zone = np.array([[1, 1], [1, 1]])  # 4-pixel zone
        >>> density = calculate_weighted_skeleton_density(skeleton, intensity, zone)
        >>> # Result: (5 + 3) / 4 = 2.0
    """
    # Extract skeleton pixels within this zone
    skeleton_in_zone = skeleton & zone_mask
    
    # Sum Sato intensity on skeleton pixels in this zone
    skeleton_intensity = np.sum(intensity_map[skeleton_in_zone])
    
    # Zone area (total pixels, not just skeleton)
    zone_area = np.sum(zone_mask)
    
    if zone_area == 0:
        return np.nan  # Can't calculate density for zero-area zone
    
    if np.sum(skeleton_in_zone) == 0:
        return 0.0  # No skeleton in this zone = zero density
    
    # Weighted density = intensity per unit area
    weighted_density = skeleton_intensity / zone_area
    
    # Sanity check
    if weighted_density < 0:
        print(f"  WARNING: Negative weighted density detected: {weighted_density}")
        return 0.0
    
    return weighted_density


def calculate_weighted_skeleton_ratios(skeleton, sato_map, zones):
    """
    Calculate ALL weighted skeleton density ratios between zones.
    
    This function computes density ratios to test the cortical fortification hypothesis:
    - High Outer/Inner ratio → Cortical enrichment (infection phenotype)
    - High Outer/Mid ratio → Peripheral vs intermediate zone comparison
    - High Mid/Inner ratio → Gradient from center outward
    
    Args:
        skeleton (ndarray): Binary skeleton (Y, X)
        sato_map (ndarray): Raw Sato tubeness (Y, X) - NOT normalized!
        zones (dict): Dict with 'inner', 'mid', 'outer' boolean masks
        
    Returns:
        dict: Contains weighted densities and all three ratios
            - 'inner_weighted_density': Zone 1 density
            - 'mid_weighted_density': Zone 2 density  
            - 'cortical_weighted_density': Zone 3 density
            - 'cortical_inner_ratio_weighted': Outer/Inner (primary metric)
            - 'cortical_mid_ratio_weighted': Outer/Mid
            - 'mid_inner_ratio_weighted': Mid/Inner
    """
    # Calculate densities for each zone
    outer_density = calculate_weighted_skeleton_density(
        skeleton, sato_map, zones['outer']
    )
    
    mid_density = calculate_weighted_skeleton_density(
        skeleton, sato_map, zones['mid']
    )
    
    inner_density = calculate_weighted_skeleton_density(
        skeleton, sato_map, zones['inner']
    )
    
    # Calculate all three ratios with proper NaN/zero handling
    
    # Ratio 1: Outer/Inner (primary cortical fortification metric)
    if np.isnan(inner_density) or inner_density == 0:
        cortical_inner_ratio = np.nan
    elif np.isnan(outer_density):
        cortical_inner_ratio = np.nan
    else:
        cortical_inner_ratio = outer_density / inner_density
    
    # Ratio 2: Outer/Mid (peripheral vs intermediate)
    if np.isnan(mid_density) or mid_density == 0:
        cortical_mid_ratio = np.nan
    elif np.isnan(outer_density):
        cortical_mid_ratio = np.nan
    else:
        cortical_mid_ratio = outer_density / mid_density
    
    # Ratio 3: Mid/Inner (intermediate vs center)
    if np.isnan(inner_density) or inner_density == 0:
        mid_inner_ratio = np.nan
    elif np.isnan(mid_density):
        mid_inner_ratio = np.nan
    else:
        mid_inner_ratio = mid_density / inner_density
    
    # Validation: Flag extreme ratios
    for ratio_name, ratio_val in [('Outer/Inner', cortical_inner_ratio), 
                                    ('Outer/Mid', cortical_mid_ratio),
                                    ('Mid/Inner', mid_inner_ratio)]:
        if not np.isnan(ratio_val):
            if ratio_val > 20 or ratio_val < 0.05:
                print(f"  WARNING: Extreme {ratio_name} ratio detected: {ratio_val:.2f}")
    
    return {
        'inner_weighted_density': inner_density,
        'mid_weighted_density': mid_density,
        'cortical_weighted_density': outer_density,
        'cortical_inner_ratio_weighted': cortical_inner_ratio,
        'cortical_mid_ratio_weighted': cortical_mid_ratio,
        'mid_inner_ratio_weighted': mid_inner_ratio
    }


# ============================================================================
# MAIN SKELETONIZER CLASS
# ============================================================================

class CellSkeletonizer:
    """Perform cell-level skeletonization with dual-path 3D processing and weighted analysis"""
    
    def __init__(self, config):
        self.full_config = config
        self.config = config['skeletonization']
        
        # V3: Dual-path processing config
        self.dual_path_config = config.get('dual_path_processing', {
            'enable': True,
            'use_3d_filtering': True,
            'path_a_reference': 'mip', # Maximum Intensity Projection
            'path_b_structure': 'sato'
        })
        
        # Parameters - UPDATED to match old pipeline defaults for better skeleton quality
        self.sato_sigmas = self.config.get('sato_sigmas', [0.7, 0.96, 1.22, 1.48, 1.74, 2.0])
        self.binary_percentile = self.config.get('binary_threshold_percentile', 55.0)
        self.min_branch_length = self.config.get('min_branch_length', 5.0)
        self.use_clahe = self.config.get('use_clahe', False)  # Disabled by default (Vieira 2025)

        # Pre-smoothing before Sato filter (old pipeline approach)
        self.pre_smoothing_sigma = self.config.get('gaussian_sigma', 0.7)

        # Advanced thresholding parameters - UPDATED to use percentile by default
        self.use_automatic_threshold = self.config.get('use_automatic_threshold', False)
        self.threshold_method = self.config.get('automatic_threshold_method', 'otsu')  # 'otsu' or 'li'
        self.use_hysteresis = self.config.get('use_hysteresis_threshold', False)
        self.hysteresis_low_fraction = self.config.get('hysteresis_low_fraction', 0.6)

        # Fidelity metric calculation
        self.calculate_fidelity = self.config.get('calculate_fidelity_metric', True)

        # Weighted skeleton analysis (ENABLED BY DEFAULT)
        self.calculate_weighted_metrics = self.config.get('calculate_weighted_skeleton_metrics', True)

        # Post-processing parameters - UPDATED to match old pipeline
        self.min_object_size = self.config.get('min_object_size', 3)
        # Minimal closing (old pipeline used 1, not aggressive 3)
        self.closing_radius = self.config.get('closing_radius', 1)
        self.opening_radius = self.config.get('opening_radius', 0)
        self.fill_holes = self.config.get('fill_holes', False)  # Should be FALSE for actin!
        self.dilation_radius = self.config.get('dilation_radius', 0)
        self.erosion_radius = self.config.get('erosion_radius', 0)
        
        # NEW: Vieira preprocessing
        self.vieira_config = config.get('vieira_preprocessing', {'enable': False, 'sigma': 0.0})
        self.vieira_sigma = self.vieira_config.get('sigma', 0.0)
        
        # Output setup
        self.save_skeletons = config['processing'].get('save_skeletons', False)
        self.save_sato_maps = config['processing'].get('save_sato_maps', False)  # NEW: For validation
        if self.save_skeletons or self.save_sato_maps:
            self.output_dir = Path(config['output_dir'])
            (self.output_dir / 'skeletons').mkdir(parents=True, exist_ok=True)
            if self.save_sato_maps:
                (self.output_dir / 'sato_maps').mkdir(parents=True, exist_ok=True)

    def _apply_sato_3d(self, volume):
        """Apply Sato filter on 3D volume (Z, Y, X)"""
        if volume.ndim != 3:
            return filters.sato(volume, sigmas=self.sato_sigmas, mode='reflect')
        return filters.sato(volume, sigmas=self.sato_sigmas, mode='reflect')
    
    def _calculate_intensity_capture_ratio(self, skeleton, intensity_image, mask):
        """
        Calculate Intensity Capture Ratio (Fidelity Metric)
        
        Higher ratios indicate better skeleton fidelity to the actual actin structure.
        """
        if np.sum(skeleton) == 0:
            return 0.0
        
        # Dilate skeleton by 1 pixel to create capture zone
        dilated_skeleton = morphology.binary_dilation(skeleton, morphology.disk(1))
        
        # Calculate intensity sums
        total_intensity = np.sum(intensity_image[mask])
        captured_intensity = np.sum(intensity_image[dilated_skeleton & mask])
        
        if total_intensity == 0:
            return 0.0
            
        return captured_intensity / total_intensity
    
    def _automatic_threshold(self, image, mask, method='otsu'):
        """Calculate automatic threshold using Otsu or Li's method"""
        valid_pixels = image[mask]
        
        if len(valid_pixels) == 0:
            return 0.0
        
        try:
            if method.lower() == 'li':
                threshold = filters.threshold_li(valid_pixels)
            else:  # Default to Otsu
                threshold = filters.threshold_otsu(valid_pixels)
            return threshold
        except Exception as e:
            print(f"  Warning: Automatic thresholding failed ({str(e)}), using percentile fallback")
            return np.percentile(valid_pixels, self.binary_percentile)
    
    def _hysteresis_threshold(self, image, mask, high_threshold, low_fraction=0.5):
        """Apply hysteresis thresholding for better structural preservation"""
        low_threshold = high_threshold * low_fraction
        binary = filters.apply_hysteresis_threshold(image, low_threshold, high_threshold)
        binary = binary & mask
        return binary

    def skeletonize_cell(self, cell_mask, actin_image, image_name, cell_id,
                         nucleus_mask=None, nuclei_image=None, actin_volume_3d=None,
                         actin_image_raw=None, actin_volume_3d_raw=None, **kwargs):
        """
        Main processing function for a single cell with weighted skeleton analysis.

        Parameters:
        -----------
        actin_image_raw : ndarray, optional
            Raw (unadjusted) 2D actin image for B/C/G-independent measurements.
            Falls back to actin_image if not provided (backward compatible).
        actin_volume_3d_raw : ndarray, optional
            Raw (unadjusted) 3D actin volume for B/C/G-independent measurements.
            Falls back to actin_volume_3d if not provided.
        """
        try:
            # 1. Extract Cell ROI (Bounding Box)
            slices = ndimage.find_objects(cell_mask)
            if not slices:
                return None
            bbox = slices[0]

            # Crop 2D data
            mask_crop = cell_mask[bbox]
            actin_crop = actin_image[bbox]

            # Raw data for measurements (fallback to preprocessed for backward compat)
            if actin_image_raw is None:
                actin_image_raw = actin_image
            if actin_volume_3d_raw is None:
                actin_volume_3d_raw = actin_volume_3d
            actin_crop_raw = actin_image_raw[bbox]

            # Crop Nucleus Mask if available
            if nucleus_mask is not None:
                nuc_mask_crop = nucleus_mask[bbox]
                nuc_mask_crop = (nuc_mask_crop > 0)
            else:
                nuc_mask_crop = None

            # Determine Processing Mode (2D vs 3D)
            use_3d = (actin_volume_3d is not None) and \
                     (self.dual_path_config['use_3d_filtering']) and \
                     (actin_volume_3d.ndim == 3) and \
                     (actin_volume_3d.shape[0] > 1)

            processing_mode = "3D" if use_3d else "2D"

            # ============================================================
            # PATH A: REFERENCE MAP (Biological Mass - RAW data)
            # ============================================================
            # Uses RAW (unadjusted) data so measurements are independent of
            # brightness/contrast/gamma settings used for Cellpose segmentation.
            if use_3d and actin_volume_3d_raw is not None:
                z_crop_raw = actin_volume_3d_raw[:, bbox[0], bbox[1]]
                intensity_ref = np.max(z_crop_raw, axis=0)  # MIP of raw
            else:
                intensity_ref = actin_crop_raw.copy()

            intensity_ref = intensity_ref * mask_crop

            # ============================================================
            # PATH B: STRUCTURE MAP (Topology Detection - Sato on RAW)
            # ============================================================
            # Sato is computed on RAW data. The per-cell normalization (below)
            # removes scale effects, so skeleton topology is identical to
            # preprocessed. But sato_raw now reflects true biological signal.

            # B0. Pre-smoothing BEFORE Sato
            if self.pre_smoothing_sigma > 0:
                actin_crop_smoothed = ndimage.gaussian_filter(actin_crop_raw.astype(float), sigma=self.pre_smoothing_sigma)
            else:
                actin_crop_smoothed = actin_crop_raw.astype(float)

            # B1. Feature Detection (Sato)
            if use_3d and actin_volume_3d_raw is not None:
                cell_vol_raw = actin_volume_3d_raw[:, bbox[0], bbox[1]]
                if self.pre_smoothing_sigma > 0:
                    cell_vol_smoothed = ndimage.gaussian_filter(cell_vol_raw.astype(float), sigma=self.pre_smoothing_sigma)
                else:
                    cell_vol_smoothed = cell_vol_raw.astype(float)
                tubeness_3d = self._apply_sato_3d(cell_vol_smoothed)
                structure_map = np.max(tubeness_3d, axis=0)
            else:
                structure_map = filters.sato(actin_crop_smoothed, sigmas=self.sato_sigmas, mode='reflect')

            # *** CRITICAL: Store raw Sato BEFORE normalization ***
            # Now derived from raw data — B/C/G-independent for weighted analysis
            sato_raw = structure_map.copy()

            # B2. Normalization (Per-Cell - ONLY for binarization)
            sm_min, sm_max = structure_map.min(), structure_map.max()
            if sm_max > sm_min:
                norm_map = (structure_map - sm_min) / (sm_max - sm_min)
            else:
                norm_map = structure_map

            # B3. CLAHE Enhancement - DISABLED by default (Vieira 2025 recommendation)
            if self.use_clahe:
                norm_map = exposure.equalize_adapthist(norm_map, clip_limit=0.01)

            # B4. Vieira Smoothing - Optional light smoothing
            if self.vieira_sigma > 0:
                norm_map = filters.gaussian(norm_map, sigma=self.vieira_sigma)
                print(f"  Cell {cell_id}: Applied Vieira smoothing (σ={self.vieira_sigma})")

            # B5. Binarization with Advanced Thresholding
            valid_pixels = norm_map[mask_crop > 0]
            if len(valid_pixels) == 0:
                return None
            
            # Determine threshold
            if self.use_automatic_threshold:
                thresh = self._automatic_threshold(norm_map, mask_crop, method=self.threshold_method)
                print(f"  Cell {cell_id}: Auto-threshold ({self.threshold_method}) = {thresh:.4f}")
            else:
                thresh = np.percentile(valid_pixels, self.binary_percentile)
                print(f"  Cell {cell_id}: Manual threshold ({self.binary_percentile}th percentile) = {thresh:.4f}")
            
            # Apply thresholding
            if self.use_hysteresis and self.use_automatic_threshold:
                binary = self._hysteresis_threshold(
                    norm_map, mask_crop, high_threshold=thresh,
                    low_fraction=self.hysteresis_low_fraction
                )
                print(f"  Cell {cell_id}: Applied hysteresis (low={thresh * self.hysteresis_low_fraction:.4f})")
            else:
                binary = norm_map > thresh
                binary = binary & mask_crop
            
            # B6. Morphological Post-Processing
            if self.min_object_size > 0:
                binary = remove_small_objects(binary, min_size=self.min_object_size)
            
            if self.closing_radius > 0:
                binary = binary_closing(binary, disk(self.closing_radius))
                
            if self.opening_radius > 0:
                binary = morphology.binary_opening(binary, disk(self.opening_radius))
            
            if self.fill_holes:
                binary = morphology.remove_small_holes(binary, area_threshold=64)
            
            if self.dilation_radius > 0:
                binary = morphology.binary_dilation(binary, disk(self.dilation_radius))
            
            if self.erosion_radius > 0:
                binary = morphology.binary_erosion(binary, disk(self.erosion_radius))
            
            binary = binary & mask_crop

            # B7. Skeletonization
            skeleton = skeletonize(binary)
            
            # ============================================================
            # GRAPH CONSTRUCTION & ANALYSIS (Weighted with Sato!)
            # ============================================================
            
            if np.sum(skeleton) < 5:
                return None

            # *** KEY CHANGE: Use sato_raw for weighting (not intensity_ref) ***
            skeleton_obj = Skeleton(skeleton, source_image=sato_raw)
            graph_summary = summarize(skeleton_obj)
            
            num_branches = len(graph_summary)
            num_junctions = max(0, num_branches - 1)
            
            # Traditional topology metrics
            features = {
                'total_length': graph_summary['branch-distance'].sum(),
                'num_branches': num_branches,
                'mean_branch_length': graph_summary['branch-distance'].mean() if num_branches > 0 else 0,
                'num_junctions': num_junctions,
                'num_endpoints': max(0, num_branches * 2 - (num_junctions * 3)),
                'skeleton_density': np.sum(skeleton) / np.sum(mask_crop)
            }
            
            # *** NEW: Extract weighted metrics from graph_summary ***
            if 'mean-pixel-value' in graph_summary.columns:
                features['total_network_intensity'] = float(graph_summary['mean-pixel-value'].sum())
                features['mean_bundle_intensity'] = float(graph_summary['mean-pixel-value'].mean())
                features['weighted_total_length'] = float((graph_summary['branch-distance'] * 
                                                           graph_summary['mean-pixel-value']).sum())
                print(f"  Cell {cell_id}: Mean bundle intensity (Sato) = {features['mean_bundle_intensity']:.2f}")
            else:
                features['total_network_intensity'] = np.nan
                features['mean_bundle_intensity'] = np.nan
                features['weighted_total_length'] = np.nan
            
            # Calculate Fidelity Metric
            if self.calculate_fidelity:
                intensity_capture_ratio = self._calculate_intensity_capture_ratio(
                    skeleton, intensity_ref, mask_crop
                )
                features['intensity_capture_ratio'] = intensity_capture_ratio
                features['fidelity'] = intensity_capture_ratio
                print(f"  Cell {cell_id}: Fidelity (Intensity Capture) = {intensity_capture_ratio:.3f}")
            else:
                features['intensity_capture_ratio'] = np.nan
                features['fidelity'] = np.nan

            # ============================================================
            # CORTICAL ANALYSIS & WEIGHTED SKELETON RATIOS
            # ============================================================
            
            cortical_results = {}
            if nuc_mask_crop is not None:
                # Legacy cortical analysis (raw intensity-based)
                eroded = morphology.binary_erosion(mask_crop, disk(5))
                cortex = mask_crop & ~eroded
                inner = mask_crop & ~cortex
                
                cortical_intensity = np.mean(intensity_ref[cortex]) if np.sum(cortex) > 0 else 0
                inner_intensity = np.mean(intensity_ref[inner]) if np.sum(inner) > 0 else 0
                
                cortical_results = {
                    'cortical_intensity': cortical_intensity,
                    'inner_intensity': inner_intensity,
                    'cortical_inner_ratio': cortical_intensity / inner_intensity if inner_intensity > 0 else 0,
                    'cortical_area': np.sum(cortex),
                    'inner_area': np.sum(inner)
                }
                
                # *** NEW: Weighted Skeleton Analysis (Sato-based) ***
                if self.calculate_weighted_metrics and TOPOLOGICAL_AVAILABLE:
                    try:
                        # Get nucleus-anchored zones (reuse topological analyzer's distance mapping)
                        topo_analyzer = TopologicalAnalyzer(self.full_config.get('topological_analysis', {}))
                        norm_map_radial, radial_method = topo_analyzer.get_normalized_distance_map(mask_crop, nuc_mask_crop)
                        
                        # Define 3 zones
                        zones = {
                            'inner': (norm_map_radial >= 0.0) & (norm_map_radial < 0.33) & mask_crop,
                            'mid':   (norm_map_radial >= 0.33) & (norm_map_radial < 0.66) & mask_crop,
                            'outer': (norm_map_radial >= 0.66) & mask_crop
                        }
                        
                        # Calculate weighted skeleton densities and ratios
                        weighted_results = calculate_weighted_skeleton_ratios(skeleton, sato_raw, zones)
                        cortical_results.update(weighted_results)
                        
                        print(f"  Cell {cell_id}: Weighted ratios - Outer/Inner={weighted_results['cortical_inner_ratio_weighted']:.2f}, "
                              f"Outer/Mid={weighted_results['cortical_mid_ratio_weighted']:.2f}, "
                              f"Mid/Inner={weighted_results['mid_inner_ratio_weighted']:.2f}")
                    
                    except Exception as e:
                        print(f"  WARNING: Weighted skeleton analysis failed for cell {cell_id}: {str(e)}")
                        # Add NaN placeholders
                        cortical_results.update({
                            'inner_weighted_density': np.nan,
                            'mid_weighted_density': np.nan,
                            'cortical_weighted_density': np.nan,
                            'cortical_inner_ratio_weighted': np.nan,
                            'cortical_mid_ratio_weighted': np.nan,
                            'mid_inner_ratio_weighted': np.nan
                        })

            # ============================================================
            # TOPOLOGICAL ANALYSIS
            # ============================================================
            
            topological_features = None
            radial_method = "none"
            
            if TOPOLOGICAL_AVAILABLE:
                topo_analyzer = TopologicalAnalyzer(self.full_config.get('topological_analysis', {}))
                
                topo_metrics = topo_analyzer.analyze_topology(
                    skeleton_data={
                        'skeleton': skeleton,
                        'graph_obj': skeleton_obj,
                        'intensity_ref': intensity_ref
                    },
                    cell_mask=mask_crop,
                    nucleus_mask=nuc_mask_crop
                )
                
                if topo_metrics:
                    topological_features = topo_metrics
                    radial_method = topo_metrics.get('radial_method', 'unknown')
                else:
                    topological_features = self._get_empty_topo_features()

            # ============================================================
            # SAVE & RETURN
            # ============================================================
            
            if self.save_skeletons:
                skeleton_dir = self.output_dir / 'skeletons'
                full_shape = cell_mask.shape
                skeleton_full = np.zeros(full_shape, dtype=np.uint8)
                skeleton_full[bbox] = skeleton.astype(np.uint8) * 255
                
                tifffile.imwrite(
                    skeleton_dir / f'cell_{cell_id:03d}_skeleton.tif',
                    skeleton_full
                )
            
            # Save Sato map for validation
            if self.save_sato_maps:
                sato_dir = self.output_dir / 'sato_maps'
                # Normalize for saving (0-255)
                sato_normalized = (sato_raw - sato_raw.min()) / (sato_raw.max() - sato_raw.min() + 1e-9)
                sato_full = np.zeros(cell_mask.shape, dtype=np.uint8)
                sato_full[bbox] = (sato_normalized * 255).astype(np.uint8)
                
                tifffile.imwrite(
                    sato_dir / f'cell_{cell_id:03d}_sato.tif',
                    sato_full
                )
            
            return {
                'skeleton': skeleton,
                'sato_raw': sato_raw,  # NEW: For downstream analysis/visualization
                'features': features,
                'cortical_features': cortical_results,
                'topological_features': topological_features,
                'radial_method': radial_method,
                'graph_obj': skeleton_obj,
                'graph_summary': graph_summary,
                'bbox': bbox,
                'processing_mode': processing_mode
            }
            
        except Exception as e:
            print(f"  ERROR: Skeletonization failed for cell {cell_id}: {str(e)}")
            import traceback
            traceback.print_exc()
            return None
    
    def _get_empty_topo_features(self):
        """Return empty topological features dict"""
        return {
            'zone_1_density': np.nan,
            'zone_2_density': np.nan,
            'zone_3_density': np.nan,
            'peripheral_gradient': np.nan,
            'spectral_entropy': np.nan,
            'zone_3_lacunarity': np.nan
        }