"""
Filtering Module for Bacterial Burden Pipeline
Hierarchical filtering: edge clearing -> cell morphology -> nucleus association -> nucleus morphology
"""

import numpy as np
from skimage import measure, segmentation


class HostCellFilter:
    """
    Filter host cells for bacterial burden analysis.

    Four-stage filtering hierarchy:
    1. clear_border: remove cells touching the image edge
    2. Cell morphology: area, eccentricity, solidity thresholds
    3. Parent-child: keep only cells with exactly 1 nucleus
    4. Nucleus morphology: area, eccentricity, solidity, nuc/cell ratio
       (catches dividing cells and bacteria-damaged nuclei)
    """

    def __init__(self, config):
        self.filter_config = config['filtering']
        self.host_config = config['filtering']['host_cell']
        self.nucleus_config = config['filtering'].get('nucleus', {})
        self.require_single_nucleus = config['filtering'].get('require_single_nucleus', True)

    def _find_nucleus_in_cell(self, cell_mask_binary, nucleus_masks):
        """
        Find nuclei overlapping with a host cell using 50% overlap rule.

        A nucleus is considered "inside" the cell if >50% of its pixels
        overlap with the cell mask.

        Parameters
        ----------
        cell_mask_binary : ndarray, bool [Y, X]
            Binary mask for a single host cell
        nucleus_masks : ndarray, int [Y, X]
            Labeled nucleus mask image

        Returns
        -------
        tuple (primary_nucleus_id, valid_count)
            primary_nucleus_id: label of the first valid nucleus, or None
            valid_count: number of nuclei meeting the overlap criterion
        """
        nuclei_in_cell = nucleus_masks * cell_mask_binary
        unique_nuclei = np.unique(nuclei_in_cell)
        unique_nuclei = unique_nuclei[unique_nuclei > 0]

        valid_count = 0
        primary_nucleus = None

        for nuc_id in unique_nuclei:
            overlap = np.sum(nuclei_in_cell == nuc_id)
            nucleus_area = np.sum(nucleus_masks == nuc_id)
            if nucleus_area > 0 and (overlap / nucleus_area) > 0.5:
                valid_count += 1
                if primary_nucleus is None:
                    primary_nucleus = nuc_id

        return primary_nucleus, valid_count

    def filter_host_cells(self, host_masks, nucleus_masks):
        """
        Apply hierarchical filtering to host cell masks.

        Parameters
        ----------
        host_masks : ndarray, int [Y, X]
            Labeled host cell masks from segmentation
        nucleus_masks : ndarray, int [Y, X]
            Labeled nucleus masks from segmentation

        Returns
        -------
        dict
            valid_cells : list of dict
                Each dict has: cell_id, nucleus_id, properties (cell_area, etc.)
            excluded_cells : list of dict
                Each dict has: cell_id, reason
            filtered_host_masks : ndarray
                Host masks with excluded cells zeroed out
            exclusion_stats : dict
                Summary counts by exclusion reason
        """
        # Stage 1: Clear border cells
        if self.filter_config.get('clear_border', True):
            host_masks = segmentation.clear_border(host_masks)

        # Get region properties
        cell_props = measure.regionprops(host_masks)

        valid_cells = []
        excluded_cells = []
        excluded_ids = set()

        for prop in cell_props:
            cell_id = prop.label
            exclude_reason = None

            # Stage 2: Morphology filters
            if prop.area < self.host_config['min_area']:
                exclude_reason = 'too_small'
            elif prop.area > self.host_config['max_area']:
                exclude_reason = 'too_large'
            elif prop.eccentricity > self.host_config.get('max_eccentricity', 0.95):
                exclude_reason = 'too_elongated'
            elif prop.solidity < self.host_config.get('min_solidity', 0.7):
                exclude_reason = 'irregular_shape'

            # Stage 3: Parent-child nucleus association
            nuc_id = None
            if exclude_reason is None:
                cell_binary = (host_masks == cell_id)
                nuc_id, nuc_count = self._find_nucleus_in_cell(cell_binary, nucleus_masks)

                if self.require_single_nucleus:
                    if nuc_id is None or nuc_count == 0:
                        exclude_reason = 'no_nucleus'
                    elif nuc_count > 1:
                        exclude_reason = 'multiple_nuclei'

            # Stage 4: Nucleus morphology filtering
            nuc_props = None
            if exclude_reason is None and nuc_id is not None and self.nucleus_config:
                nuc_binary = (nucleus_masks == nuc_id)
                nuc_regions = measure.regionprops(nuc_binary.astype(int))

                if nuc_regions:
                    nuc_props = nuc_regions[0]
                    nuc_area = nuc_props.area
                    nuc_cell_ratio = nuc_area / prop.area if prop.area > 0 else 0.0

                    if nuc_area < self.nucleus_config.get('min_area', 0):
                        exclude_reason = 'nucleus_too_small'
                    elif nuc_area > self.nucleus_config.get('max_area', float('inf')):
                        exclude_reason = 'nucleus_too_large'
                    elif nuc_props.eccentricity > self.nucleus_config.get('max_eccentricity', 1.0):
                        exclude_reason = 'nucleus_too_elongated'
                    elif nuc_props.solidity < self.nucleus_config.get('min_solidity', 0.0):
                        exclude_reason = 'nucleus_irregular_shape'
                    elif nuc_cell_ratio < self.nucleus_config.get('min_nuc_cell_ratio', 0.0):
                        exclude_reason = 'nucleus_ratio_abnormal'
                    elif nuc_cell_ratio > self.nucleus_config.get('max_nuc_cell_ratio', 1.0):
                        exclude_reason = 'nucleus_ratio_abnormal'

            if exclude_reason is not None:
                excluded_cells.append({
                    'cell_id': cell_id,
                    'reason': exclude_reason,
                    'area': prop.area,
                    'eccentricity': prop.eccentricity,
                    'solidity': prop.solidity
                })
                excluded_ids.add(cell_id)
            else:
                # Build properties dict with nucleus info
                cell_properties = {
                    'cell_area': prop.area,
                    'cell_solidity': prop.solidity,
                    'cell_eccentricity': prop.eccentricity,
                    'cell_perimeter': prop.perimeter
                }
                if nuc_props is not None:
                    cell_properties['nucleus_area'] = nuc_props.area
                    cell_properties['nucleus_solidity'] = nuc_props.solidity
                    cell_properties['nucleus_eccentricity'] = nuc_props.eccentricity
                    cell_properties['nuc_cell_ratio'] = nuc_props.area / prop.area if prop.area > 0 else 0.0

                valid_cells.append({
                    'cell_id': cell_id,
                    'nucleus_id': nuc_id,
                    'properties': cell_properties
                })

        # Create filtered mask
        filtered_masks = host_masks.copy()
        for cid in excluded_ids:
            filtered_masks[host_masks == cid] = 0

        # Exclusion statistics
        reason_counts = {}
        for exc in excluded_cells:
            r = exc['reason']
            reason_counts[r] = reason_counts.get(r, 0) + 1

        total = len(cell_props)
        exclusion_stats = {
            'total_detected': total,
            'valid_count': len(valid_cells),
            'excluded_count': len(excluded_cells),
            'exclusion_rate': len(excluded_cells) / total if total > 0 else 0.0,
            'reasons': reason_counts
        }

        return {
            'valid_cells': valid_cells,
            'excluded_cells': excluded_cells,
            'filtered_host_masks': filtered_masks,
            'exclusion_stats': exclusion_stats
        }
