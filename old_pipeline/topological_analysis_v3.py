"""
Topological Analysis Module V3 - 3-Zone Fortification Metrics
"""

import numpy as np
from scipy import ndimage
import warnings

class TopologicalAnalyzer:
    """
    Quantifies "Cortical Fortification" using a 3-Zone Nucleus-Anchored Model.
    """
    
    def __init__(self, config=None):
        self.config = config or {}
        
        # Handle config format safely
        if isinstance(self.config, dict):
            self.n_zones = self.config.get('n_zones', 3)
        else:
            self.n_zones = 3
            
        # Zone definitions (Normalized Distance 0.0=Nucleus -> 1.0=Membrane)
        self.zone_thresholds = {
            'inner': 0.33,
            'mid': 0.66
        }

    def get_normalized_distance_map(self, cell_mask, nucleus_mask):
        """
        Creates a relative distance map (0.0 at nucleus edge, 1.0 at cell edge).
        Handles irregular shapes and off-center nuclei.
        """
        # Safety: If no nucleus, fallback to centroid distance (approx)
        if nucleus_mask is None or np.sum(nucleus_mask) == 0:
            dist_to_mem = ndimage.distance_transform_edt(cell_mask)
            max_dist = np.max(dist_to_mem)
            if max_dist == 0: return np.zeros_like(dist_to_mem), "error"
            # Inverted: max at center (0), 0 at edge
            # We want 0 at center, 1 at edge
            return 1.0 - (dist_to_mem / (max_dist + 1e-9)), "centroid_fallback"

        # 1. Distance FROM Nucleus (moving outward)
        # Invert nucleus: 0 inside nucleus, 1 outside
        # Transform calculates distance to nearest 0 (nucleus)
        dist_from_nuc = ndimage.distance_transform_edt(1 - nucleus_mask)
        
        # 2. Distance FROM Membrane (moving inward)
        # Transform calculates distance to nearest 0 (background)
        dist_from_mem = ndimage.distance_transform_edt(cell_mask)

        # 3. Normalized Ratio
        # Map = D_nuc / (D_nuc + D_mem)
        # At nucleus edge: D_nuc=0 -> Map=0
        # At membrane: D_mem=0 -> Map=1
        total_dist = dist_from_nuc + dist_from_mem
        
        # Avoid division by zero
        with np.errstate(divide='ignore', invalid='ignore'):
            normalized_map = np.divide(dist_from_nuc, total_dist, out=np.zeros_like(total_dist), where=total_dist!=0)
        
        # Mask outside cell (ensure background is NaN or handled)
        normalized_map = np.where(cell_mask, normalized_map, np.nan)
        
        return normalized_map, "nucleus_anchored"

    def analyze_topology(self, skeleton_data, cell_mask, nucleus_mask=None):
        """
        Calculates 3-zone metrics: Density, Lacunarity, Gradient, Entropy.
        
        Args:
            skeleton_data (dict): MUST contain 'graph_obj', 'intensity_ref', 'skeleton'
            cell_mask (ndarray): Boolean mask of cell
            nucleus_mask (ndarray): Boolean mask of nucleus
        """
        if skeleton_data is None: 
            return None

        # Unpack data
        skel_obj = skeleton_data.get('graph_obj')
        intensity_ref = skeleton_data.get('intensity_ref')
        skeleton_pixels = skeleton_data.get('skeleton') # Binary
        
        if skel_obj is None:
            return None
        
        # 1. Generate Map
        norm_map, method = self.get_normalized_distance_map(cell_mask, nucleus_mask)
        
        # 2. Define Zones
        # Zone 1: Inner (< 0.33)
        # Zone 2: Mid (0.33 - 0.66)
        # Zone 3: Outer (> 0.66)
        
        # Handle NaNs in map
        safe_map = np.nan_to_num(norm_map, nan=-1.0)
        
        zones = {
            'inner': (safe_map >= 0.0) & (safe_map < self.zone_thresholds['inner']),
            'mid':   (safe_map >= self.zone_thresholds['inner']) & (safe_map < self.zone_thresholds['mid']),
            'outer': (safe_map >= self.zone_thresholds['mid'])
        }
        
        # Enforce cell mask on zones
        for k in zones:
            zones[k] = zones[k] & cell_mask

        # 3. Calculate Metrics
        metrics = {'radial_method': method}
        
        # A. Node Density (Junctions per pixel)
        # skan 0.10+ uses degrees and coordinates
        if hasattr(skel_obj, 'degrees') and hasattr(skel_obj, 'coordinates'):
            degrees = skel_obj.degrees
            coords = skel_obj.coordinates.astype(int)
            
            node_counts = {'inner': 0, 'mid': 0, 'outer': 0}
            
            # Iterate nodes
            for i, deg in enumerate(degrees):
                if deg > 2: # Junction
                    # Coords are (row, col) usually
                    if coords.shape[1] >= 2:
                        r, c = coords[i, 0], coords[i, 1]
                        
                        # Check bounds
                        if 0 <= r < norm_map.shape[0] and 0 <= c < norm_map.shape[1]:
                            val = safe_map[r, c]
                            if val == -1: continue
                            
                            if val < self.zone_thresholds['inner']: node_counts['inner'] += 1
                            elif val < self.zone_thresholds['mid']: node_counts['mid'] += 1
                            else: node_counts['outer'] += 1
        else:
            # Fallback if skan version differs or graph empty
            node_counts = {'inner': 0, 'mid': 0, 'outer': 0}

        # Calculate density
        for z_name, mask in zones.items():
            area = np.sum(mask)
            if area > 0:
                metrics[f'zone_{z_name[-1]}_density'] = node_counts[z_name] / area # For backward compat keys
                metrics[f'{z_name}_density'] = node_counts[z_name] / area
            else:
                metrics[f'zone_{z_name[-1]}_density'] = np.nan
                metrics[f'{z_name}_density'] = np.nan

        # Map to old keys expected by CSV: zone_1 (inner), zone_2 (mid), zone_3 (outer)
        metrics['zone_1_density'] = metrics['inner_density']
        metrics['zone_2_density'] = metrics['mid_density']
        metrics['zone_3_density'] = metrics['outer_density']

        # B. Robust Peripheral Gradient
        # (Outer - Mid) -> Positive = Fortification
        if not np.isnan(metrics['outer_density']) and not np.isnan(metrics['mid_density']):
            metrics['peripheral_gradient'] = metrics['outer_density'] - metrics['mid_density']
        else:
            metrics['peripheral_gradient'] = np.nan

        # C. Zonal Lacunarity (Hole Size)
        # Invert Skeleton
        inverted_skel = (~skeleton_pixels) & cell_mask
        dist_transform = ndimage.distance_transform_edt(inverted_skel)
        
        # Mean hole size in Outer vs Inner
        for z_key, z_mask in [('outer', zones['outer']), ('inner', zones['inner'])]:
            if np.sum(z_mask) > 0:
                metrics[f'zone_{"3" if z_key=="outer" else "1"}_lacunarity'] = np.mean(dist_transform[z_mask])
            else:
                metrics[f'zone_{"3" if z_key=="outer" else "1"}_lacunarity'] = np.nan

        # D. Spectral Entropy (Chaos) in Cortex
        # Extract pixels from intensity_ref in Outer Zone
        cortex_pixels = intensity_ref[zones['outer']]
        if cortex_pixels.size > 50:
            try:
                # 1D FFT on flattened cortex pixels
                fft_vals = np.fft.fft(cortex_pixels)
                power_spectrum = np.abs(fft_vals)**2
                # Normalize to probability distribution
                ps_sum = np.sum(power_spectrum)
                if ps_sum > 0:
                    ps_norm = power_spectrum / ps_sum
                    # Shannon Entropy
                    entropy = -np.sum(ps_norm * np.log2(ps_norm + 1e-12))
                    metrics['spectral_entropy'] = entropy
                else:
                    metrics['spectral_entropy'] = 0
            except:
                metrics['spectral_entropy'] = np.nan
        else:
            metrics['spectral_entropy'] = np.nan

        return metrics