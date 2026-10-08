"""
Filtering Module v2 - Enhanced with Phenotype Rescue & Quality Tracking
Addresses the "Filtering Paradox" - prevents removal of strongest phenotypes

Changes from v1:
- Added phenotype vs quality exclusion tracking
- Added rescue mode for relaxed shape-based filtering
- Added detailed exclusion statistics
- Compatible with original config structure
"""

import numpy as np
from skimage import measure
from scipy import ndimage


class NucleiFilter:
    """
    Filter nuclei based on quality metrics to exclude irregular nuclei.
    
    NEW in v2: Distinguishes between:
    - Quality exclusions (artifacts, debris) - always remove
    - Phenotype exclusions (shape changes from biology) - optionally rescue
    """
    
    # Define which exclusion reasons are phenotype-related vs quality-related
    PHENOTYPE_EXCLUSIONS = {'irregular_shape', 'too_elongated', 'irregular_boundary'}
    QUALITY_EXCLUSIONS = {'too_small', 'too_large', 'low_intensity', 'high_variance_bacteria', 'too_close'}
    
    def __init__(self, config):
        self.config = config['nucleus_filtering']
        
        # Rescue mode settings (new in v2)
        self.rescue_mode = self.config.get('rescue_phenotype_exclusions', False)
        self.rescue_solidity = self.config.get('rescue_solidity_threshold', 0.5)
        self.rescue_eccentricity = self.config.get('rescue_eccentricity_threshold', 0.98)
    
    def calculate_roundness(self, region):
        """
        Calculate roundness: 4Ï€*Area / PerimeterÂ²
        Perfect circle = 1.0, lower values = more elongated
        """
        if region.perimeter == 0:
            return 0
        return (4 * np.pi * region.area) / (region.perimeter ** 2)
    
    def calculate_intensity_metrics(self, region, image):
        """Calculate mean intensity and std deviation for region"""
        coords = region.coords
        intensities = image[coords[:, 0], coords[:, 1]]
        
        # Normalize to 0-1
        intensities = intensities.astype(np.float32)
        if intensities.max() > 0:
            intensities = intensities / intensities.max()
        
        return {
            'mean_intensity': np.mean(intensities),
            'std_intensity': np.std(intensities)
        }
    
    def _get_effective_thresholds(self):
        """
        Get effective thresholds based on rescue mode.
        Returns relaxed thresholds if rescue mode is enabled.
        """
        if self.rescue_mode:
            return {
                'min_solidity': self.rescue_solidity,
                'max_eccentricity': self.rescue_eccentricity,
                'min_roundness': self.config.get('min_roundness', 0.4) * 0.5  # Relax by 50%
            }
        else:
            return {
                'min_solidity': self.config['min_solidity'],
                'max_eccentricity': self.config['max_eccentricity'],
                'min_roundness': self.config['min_roundness']
            }
    
    def filter_nuclei(self, nucleus_masks, nuclei_image, image_name):
        """
        Filter nuclei based on quality metrics.
        
        Parameters:
        -----------
        nucleus_masks : ndarray
            Nucleus segmentation masks
        nuclei_image : ndarray
            Raw nuclei image for intensity analysis
        image_name : str
            Name of source image (for logging)
        
        Returns:
        --------
        result : dict
            Dictionary containing:
            - valid_nuclei: List of valid nucleus IDs
            - excluded_nuclei: List of dicts with exclusion reasons and types
            - filtered_masks: Masks with only valid nuclei
            - exclusion_stats: Summary statistics (NEW in v2)
            - rescued_nuclei: List of nuclei that would have been excluded but were rescued (NEW)
        """
        nucleus_props = measure.regionprops(nucleus_masks, intensity_image=nuclei_image)
        
        valid_nuclei = []
        excluded_nuclei = []
        rescued_nuclei = []  # NEW: Track rescued phenotype exclusions
        
        # Get effective thresholds (may be relaxed in rescue mode)
        thresholds = self._get_effective_thresholds()
        
        for prop in nucleus_props:
            nucleus_id = prop.label
            exclude_reason = None
            exclusion_type = None
            would_rescue = False
            
            # 1. Check area (quality criterion - never rescue)
            if prop.area < self.config['min_area']:
                exclude_reason = 'too_small'
                exclusion_type = 'quality'
            elif prop.area > self.config['max_area']:
                exclude_reason = 'too_large'
                exclusion_type = 'quality'
            
            # 2. Check solidity (phenotype criterion - can rescue)
            elif prop.solidity < thresholds['min_solidity']:
                exclude_reason = 'irregular_shape'
                exclusion_type = 'phenotype'
                # Check if this would have been excluded without rescue
                if self.rescue_mode and prop.solidity >= self.rescue_solidity:
                    would_rescue = True
            
            # 3. Check eccentricity
            elif prop.eccentricity < self.config.get('min_eccentricity', 0.0):
                exclude_reason = 'too_round'
                exclusion_type = 'phenotype'
            elif prop.eccentricity > thresholds['max_eccentricity']:
                exclude_reason = 'too_elongated'
                exclusion_type = 'phenotype'
                if self.rescue_mode and prop.eccentricity <= self.rescue_eccentricity:
                    would_rescue = True
            
            # 4. Check roundness (phenotype criterion)
            else:
                roundness = self.calculate_roundness(prop)
                if roundness < thresholds['min_roundness']:
                    exclude_reason = 'irregular_boundary'
                    exclusion_type = 'phenotype'
                    if self.rescue_mode:
                        would_rescue = True
            
            # 5. Check intensity (quality criterion - bacteria detection)
            if exclude_reason is None:
                intensity_metrics = self.calculate_intensity_metrics(prop, nuclei_image)
                
                if intensity_metrics['mean_intensity'] < self.config['min_mean_intensity']:
                    exclude_reason = 'low_intensity'
                    exclusion_type = 'quality'
                elif intensity_metrics['std_intensity'] > self.config['max_intensity_std']:
                    exclude_reason = 'high_variance_bacteria'
                    exclusion_type = 'quality'
            
            # Store result
            if exclude_reason is None:
                valid_nuclei.append(nucleus_id)
            else:
                exclusion_record = {
                    'nucleus_id': nucleus_id,
                    'reason': exclude_reason,
                    'exclusion_type': exclusion_type,  # NEW: Track type
                    'area': prop.area,
                    'solidity': prop.solidity,
                    'eccentricity': prop.eccentricity,
                    'roundness': self.calculate_roundness(prop),
                    'rescued': would_rescue  # NEW: Flag if this was rescued
                }
                
                if would_rescue:
                    # In rescue mode, keep this nucleus but flag it
                    valid_nuclei.append(nucleus_id)
                    rescued_nuclei.append(exclusion_record)
                else:
                    excluded_nuclei.append(exclusion_record)
        
        # Create filtered masks (exclude only non-rescued nuclei)
        filtered_masks = nucleus_masks.copy()
        for excluded in excluded_nuclei:
            if not excluded.get('rescued', False):
                filtered_masks[nucleus_masks == excluded['nucleus_id']] = 0

        # Proximity filter: dilate nuclei; if two overlap, exclude both
        too_close_ids = self._filter_proximity_nuclei(nucleus_masks, valid_nuclei)
        if too_close_ids:
            n_close = len(too_close_ids)
            print(f'  Proximity filter: excluded {n_close} nuclei (too_close)')
            for nuc_id in too_close_ids:
                valid_nuclei.remove(nuc_id)
                excluded_nuclei.append({
                    'nucleus_id': nuc_id,
                    'reason': 'too_close',
                    'exclusion_type': 'quality',
                    'area': int(np.sum(nucleus_masks == nuc_id)),
                    'solidity': np.nan,
                    'eccentricity': np.nan,
                    'roundness': np.nan,
                    'rescued': False
                })
                filtered_masks[nucleus_masks == nuc_id] = 0

        # Compute exclusion statistics (NEW in v2)
        exclusion_stats = self._compute_exclusion_stats(
            excluded_nuclei, rescued_nuclei, len(nucleus_props)
        )

        return {
            'valid_nuclei': valid_nuclei,
            'excluded_nuclei': excluded_nuclei,
            'filtered_masks': filtered_masks,
            'exclusion_stats': exclusion_stats,  # NEW
            'rescued_nuclei': rescued_nuclei     # NEW
        }
    
    def _compute_exclusion_stats(self, excluded, rescued, total):
        """
        Compute summary statistics for exclusions.
        Useful for detecting if filtering is too aggressive.
        """
        if total == 0:
            return {'exclusion_rate': 0, 'phenotype_rate': 0, 'rescue_rate': 0}
        
        n_excluded = len(excluded)
        n_rescued = len(rescued)
        
        # Count by type
        phenotype_excluded = sum(1 for e in excluded if e['exclusion_type'] == 'phenotype')
        quality_excluded = sum(1 for e in excluded if e['exclusion_type'] == 'quality')
        
        # Count by reason
        reason_counts = {}
        for e in excluded + rescued:
            reason = e['reason']
            reason_counts[reason] = reason_counts.get(reason, 0) + 1
        
        return {
            'total_nuclei': total,
            'excluded_count': n_excluded,
            'rescued_count': n_rescued,
            'exclusion_rate': n_excluded / total,
            'phenotype_exclusion_rate': phenotype_excluded / total,
            'quality_exclusion_rate': quality_excluded / total,
            'rescue_rate': n_rescued / (phenotype_excluded + n_rescued) if (phenotype_excluded + n_rescued) > 0 else 0,
            'reasons': reason_counts
        }

    def _filter_proximity_nuclei(self, nucleus_masks, valid_nuclei):
        """
        Dilate each valid nucleus by proximity_dilation_pixels.
        If any two dilated masks overlap, all involved nuclei are flagged as 'too_close'.

        Uses an O(N) count map: sum dilated masks pixel-wise; pixels >= 2 mean overlap.
        Efficient for dense fields — no pairwise comparison needed.

        Returns
        -------
        set of nucleus_ids to exclude due to proximity
        """
        prox_cfg = self.config.get('proximity_exclusion', {})
        if not prox_cfg.get('enabled', False) or len(valid_nuclei) < 2:
            return set()

        dilation_px = int(prox_cfg.get('dilation_pixels', 10))
        if dilation_px <= 0:
            return set()

        # Accumulate dilated masks into count map
        overlap_count = np.zeros(nucleus_masks.shape, dtype=np.uint8)
        dilated_per_nuc = {}
        for nuc_id in valid_nuclei:
            mask = nucleus_masks == nuc_id
            dilated = ndimage.binary_dilation(mask, iterations=dilation_px)
            dilated_per_nuc[nuc_id] = dilated
            overlap_count += dilated.astype(np.uint8)

        # Pixels touched by 2+ dilated nuclei
        overlap_pixels = overlap_count >= 2

        if not np.any(overlap_pixels):
            return set()

        # Any nucleus whose dilated mask touches an overlap pixel is too close
        too_close = set()
        for nuc_id, dilated in dilated_per_nuc.items():
            if np.any(dilated & overlap_pixels):
                too_close.add(nuc_id)

        return too_close


class CellFilter:
    """
    Filter cells based on border position and association with valid nuclei.
    Enhanced with phenotype tracking (v2).
    """
    
    PHENOTYPE_EXCLUSIONS = {'irregular_shape', 'too_elongated'}
    QUALITY_EXCLUSIONS = {'too_small', 'too_large', 'on_border', 'no_valid_nucleus', 'multiple_nuclei'}
    
    def __init__(self, config):
        self.config = config['cell_filtering']
        self.rescue_mode = self.config.get('rescue_phenotype_exclusions', False)
    
    def find_nucleus_in_cell(self, cell_mask, nucleus_masks, valid_nuclei):
        """
        Find which nucleus/nuclei overlap with this cell.
        Returns (nucleus_id, num_nuclei) or (None, 0)
        """
        cell_region = cell_mask > 0
        nuclei_in_cell = nucleus_masks * cell_region
        
        unique_nuclei = np.unique(nuclei_in_cell)
        unique_nuclei = unique_nuclei[unique_nuclei > 0]
        
        valid_count = 0
        primary_nucleus = None
        
        for nuc_id in unique_nuclei:
            if nuc_id in valid_nuclei:
                overlap = np.sum(nuclei_in_cell == nuc_id)
                nucleus_area = np.sum(nucleus_masks == nuc_id)
                
                if overlap / nucleus_area > 0.5:
                    valid_count += 1
                    if primary_nucleus is None:
                        primary_nucleus = nuc_id
        
        return primary_nucleus, valid_count
    
    def filter_cells(self, cell_masks, nucleus_masks, valid_nuclei, image_name):
        """
        Filter cells with phenotype tracking.
        
        Returns dict with:
        - valid_cells: List of dicts with cell info
        - excluded_cells: List of dicts with exclusion reasons and types
        - exclusion_stats: Summary statistics
        """
        cell_props = measure.regionprops(cell_masks)
        
        img_h, img_w = cell_masks.shape
        border_buffer = self.config['border_buffer']
        
        valid_cells = []
        excluded_cells = []
        
        for prop in cell_props:
            cell_id = prop.label
            exclude_reason = None
            exclusion_type = None
            
            # 1. Check area
            if prop.area < self.config['min_area']:
                exclude_reason = 'too_small'
                exclusion_type = 'quality'
            elif prop.area > self.config['max_area']:
                exclude_reason = 'too_large'
                exclusion_type = 'quality'
            
            # 2. Check solidity (phenotype)
            elif prop.solidity < self.config['min_solidity']:
                exclude_reason = 'irregular_shape'
                exclusion_type = 'phenotype'
            
            # 3. Check eccentricity (phenotype)
            elif 'max_eccentricity' in self.config and prop.eccentricity > self.config['max_eccentricity']:
                exclude_reason = 'too_elongated'
                exclusion_type = 'phenotype'
            
            # 4. Check border (quality - never rescue)
            elif exclude_reason is None:
                bbox = prop.bbox
                if (bbox[0] < border_buffer or 
                    bbox[1] < border_buffer or
                    bbox[2] > img_h - border_buffer or 
                    bbox[3] > img_w - border_buffer):
                    exclude_reason = 'on_border'
                    exclusion_type = 'quality'
            
            # 5. Check for valid nucleus (quality)
            if exclude_reason is None:
                cell_mask = (cell_masks == cell_id)
                nucleus_id, num_nuclei = self.find_nucleus_in_cell(
                    cell_mask, nucleus_masks, valid_nuclei
                )
                
                if nucleus_id is None:
                    exclude_reason = 'no_valid_nucleus'
                    exclusion_type = 'quality'
                elif self.config.get('reject_multiple_nuclei', False) and num_nuclei > 1:
                    exclude_reason = 'multiple_nuclei'
                    exclusion_type = 'quality'
                else:
                    prop.nucleus_id = nucleus_id
                    prop.num_nuclei = num_nuclei
            
            # Store result
            if exclude_reason is None:
                valid_cells.append({
                    'cell_id': cell_id,
                    'nucleus_id': getattr(prop, 'nucleus_id', None),
                    'properties': {
                        'cell_area': prop.area,
                        'cell_solidity': prop.solidity,
                        'cell_eccentricity': prop.eccentricity,
                        'cell_perimeter': prop.perimeter
                    }
                })
            else:
                excluded_cells.append({
                    'cell_id': cell_id,
                    'reason': exclude_reason,
                    'exclusion_type': exclusion_type,
                    'area': prop.area,
                    'solidity': prop.solidity,
                    'eccentricity': prop.eccentricity
                })
        
        # Compute stats
        exclusion_stats = {
            'total_cells': len(cell_props),
            'valid_count': len(valid_cells),
            'excluded_count': len(excluded_cells),
            'exclusion_rate': len(excluded_cells) / len(cell_props) if len(cell_props) > 0 else 0
        }
        
        return {
            'valid_cells': valid_cells,
            'excluded_cells': excluded_cells,
            'exclusion_stats': exclusion_stats
        }
