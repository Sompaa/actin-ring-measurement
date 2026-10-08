"""
Bacterial Burden Quantification Module
Computes intracellular bacterial area metrics per host cell using
pre-computed Sato filter results.
"""

import numpy as np
from skimage.measure import regionprops


class BacterialBurdenQuantifier:
    """
    Quantifies bacterial burden per host cell from Sato filter area results.

    Metrics per cell:
    - host_area: pixel area of the host cell
    - intracellular_bact_area: total bacterial pixel area inside the cell
    - area_ratio: intracellular_bact_area / host_area
    """

    def __init__(self, config):
        pass

    def quantify_image(self, host_masks, valid_cells, bacteria_area_results):
        """
        Quantify all valid cells in an image.

        Parameters
        ----------
        host_masks : ndarray [Y, X]
            Filtered labeled host mask (excluded cells already zeroed)
        valid_cells : list of dict
            From HostCellFilter.filter_host_cells().
            Each dict has 'cell_id', 'nucleus_id', 'properties'.
        bacteria_area_results : dict
            Pre-computed bacterial area results from Sato filter
            (_quantify_bacteria_per_cell). Keyed by cell_id.

        Returns
        -------
        list of dict
            Per-cell measurement records with cell properties and bacterial metrics
        """
        results = []

        for cell_info in valid_cells:
            cell_id = cell_info['cell_id']

            if cell_id not in bacteria_area_results:
                continue

            area_data = bacteria_area_results[cell_id]

            cell_mask = (host_masks == cell_id)
            host_props = regionprops(cell_mask.astype(int))

            if len(host_props) == 0:
                continue

            props = host_props[0]

            # Get nucleus properties from filter output
            cell_props = cell_info.get('properties', {})

            record = {
                'cell_id': cell_id,
                'nucleus_id': cell_info.get('nucleus_id'),
                'host_area': area_data['host_area'],
                'cell_area': props.area,
                'cell_solidity': props.solidity,
                'cell_eccentricity': props.eccentricity,
                'cell_perimeter': props.perimeter,
                'nucleus_area': cell_props.get('nucleus_area'),
                'nucleus_solidity': cell_props.get('nucleus_solidity'),
                'nucleus_eccentricity': cell_props.get('nucleus_eccentricity'),
                'nuc_cell_ratio': cell_props.get('nuc_cell_ratio'),
                'intracellular_bact_area': area_data['bacterial_area'],
                'area_ratio': round(area_data['area_ratio'], 6),
            }

            results.append(record)

        return results
