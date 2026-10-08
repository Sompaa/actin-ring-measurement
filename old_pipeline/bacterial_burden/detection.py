"""
Standalone bacterial detection functions for use by both the main TopoActin pipeline
and the standalone bacterial burden module.

Extracted from segmentation_bacterial.py to allow shared usage without coupling
to the BacterialSegmenter class or its config schema.

Functions:
    normalize_cell_region       - Per-cell intensity normalization
    detect_bacteria_in_cell     - Sato ridge filter on bacteria ROI
    quantify_bacteria_per_cell  - Orchestrator: loop over cells, detect, compute area
    compute_zone_bacterial_density - Per-zone (perinuclear/middle/cortical) bacterial metrics
"""

import numpy as np


def normalize_cell_region(region, method='percentile', percentile_low=2, percentile_high=98):
    """
    Normalize intensity for a single cell's bacteria region.

    Per-cell normalization handles heterogeneous intensity distributions
    across different cells in the same image and across experiments.

    Parameters
    ----------
    region : ndarray
        2D numpy array (bacteria channel crop for one cell)
    method : str
        'percentile' or 'robust_zscore'
    percentile_low : float
        Lower percentile cutoff (default: 2)
    percentile_high : float
        Upper percentile cutoff (default: 98)

    Returns
    -------
    ndarray
        Normalized uint16 image
    """
    if region.size == 0 or region.max() == 0:
        return region

    if method == 'percentile':
        p_low = np.percentile(region, percentile_low)
        p_high = np.percentile(region, percentile_high)

        if p_high <= p_low:
            return region

        normalized = np.clip(region, p_low, p_high)
        normalized = ((normalized - p_low) / (p_high - p_low) * 65535).astype(np.uint16)

    elif method == 'robust_zscore':
        median = np.median(region)
        mad = np.median(np.abs(region - median))

        if mad == 0:
            return region

        z_score = (region.astype(np.float32) - median) / (1.4826 * mad)
        z_clipped = np.clip(z_score, -3, 3)
        normalized = ((z_clipped + 3) / 6 * 65535).astype(np.uint16)

    else:
        normalized = region

    return normalized


def detect_bacteria_in_cell(bacteria_region, params):
    """
    Detect bacterial regions using Sato ridge filter + thresholding.

    Sato filter enhances tubular/rod-shaped structures (bacteria).
    Applied to the BACTERIA channel, not actin.

    Parameters
    ----------
    bacteria_region : ndarray
        Normalized 2D array (uint16) - bacteria channel ROI
    params : dict
        Keys:
        - sato_sigmas: list of scales (e.g. [1, 2, 3])
        - sato_black_ridges: bool (False for bright bacteria on dark bg)
        - threshold_method: 'otsu', 'adaptive', or 'percentile'
        - threshold_percentile: value for 'percentile' method
        - adaptive_block_size: block size for 'adaptive' method
        - min_area: minimum bacterial region size (pixels)

    Returns
    -------
    ndarray
        Binary mask of bacterial regions
    """
    from skimage.filters import sato, threshold_otsu, threshold_local
    from skimage.morphology import remove_small_objects
    from scipy.ndimage import binary_fill_holes

    img_float = bacteria_region.astype(np.float32) / 65535.0

    sigmas = params.get('sato_sigmas', [1, 2, 3])
    black_ridges = params.get('sato_black_ridges', False)

    sato_enhanced = sato(
        img_float,
        sigmas=sigmas,
        black_ridges=black_ridges,
        mode='reflect'
    )

    threshold_method = params.get('threshold_method', 'otsu')

    if threshold_method == 'otsu':
        if sato_enhanced.max() > 0:
            thresh = threshold_otsu(sato_enhanced)
            bacteria_mask = sato_enhanced > thresh
        else:
            bacteria_mask = np.zeros_like(sato_enhanced, dtype=bool)

    elif threshold_method == 'percentile':
        percentile = params.get('threshold_percentile', 90)
        thresh = np.percentile(sato_enhanced, percentile)
        bacteria_mask = sato_enhanced > thresh

    elif threshold_method == 'adaptive':
        block_size = params.get('adaptive_block_size', 35)
        if block_size % 2 == 0:
            block_size += 1
        thresh = threshold_local(sato_enhanced, block_size=block_size, method='gaussian')
        bacteria_mask = sato_enhanced > thresh

    else:
        raise ValueError(f"Unknown threshold method: {threshold_method}")

    bacteria_mask = binary_fill_holes(bacteria_mask)

    min_area = params.get('min_area', 5)
    if min_area > 0:
        bacteria_mask = remove_small_objects(bacteria_mask, min_size=min_area)

    return bacteria_mask


def quantify_bacteria_per_cell(bacteria_2d, bacteria_raw, cell_masks, valid_cell_ids,
                                detection_params, norm_params=None):
    """
    Quantify bacterial area per cell using Sato filter + thresholding.

    Loops over all valid cells, extracts ROI, normalizes, detects bacteria
    via Sato ridge filter, and computes total bacterial pixel area.

    Parameters
    ----------
    bacteria_2d : ndarray [Y, X]
        Brightness/contrast-adjusted bacteria channel (2D MIP)
    bacteria_raw : ndarray [Y, X]
        Raw bacteria channel (before B/C) for signal range validation
    cell_masks : ndarray [Y, X]
        Labeled cell mask array (integer labels, 0=background)
    valid_cell_ids : list of int
        Cell IDs to quantify (from filtering step)
    detection_params : dict
        Bacteria detection parameters (sato_sigmas, threshold_method, min_area, etc.)
    norm_params : dict or None
        Per-cell normalization parameters. Keys:
        - enabled (bool): whether to normalize
        - method (str): 'percentile' or 'robust_zscore'
        - percentile_low (float): lower cutoff
        - percentile_high (float): upper cutoff

    Returns
    -------
    dict
        Mapping cell_id -> {
            'bacterial_area': int (total pixels),
            'host_area': int (cell pixels),
            'area_ratio': float (bacterial/host)
        }
    """
    from skimage.measure import regionprops

    if norm_params is None:
        norm_params = {
            'enabled': True,
            'method': 'percentile',
            'percentile_low': 2,
            'percentile_high': 98
        }

    padding = detection_params.get('padding', 15)
    min_range = detection_params.get('min_signal_range', 100)
    results = {}
    skipped_low_signal = 0

    print(f"Quantifying bacterial area for {len(valid_cell_ids)} cells (min_signal_range={min_range})")

    for cell_id in valid_cell_ids:
        cell_mask = (cell_masks == cell_id)
        props = regionprops(cell_mask.astype(int))

        if len(props) == 0:
            continue

        bbox = props[0].bbox
        host_area = props[0].area

        # Padded ROI
        r0 = max(0, bbox[0] - padding)
        c0 = max(0, bbox[1] - padding)
        r1 = min(bacteria_2d.shape[0], bbox[2] + padding)
        c1 = min(bacteria_2d.shape[1], bbox[3] + padding)

        cell_mask_roi = cell_mask[r0:r1, c0:c1]

        # Signal range check on RAW image (avoids false positives from B/C amplification)
        raw_check = bacteria_raw if bacteria_raw is not None else bacteria_2d
        raw_roi = raw_check[r0:r1, c0:c1]
        cell_bacteria_raw = raw_roi[cell_mask_roi]

        if cell_bacteria_raw.size > 0 and min_range > 0:
            p2 = float(np.percentile(cell_bacteria_raw, 2))
            p98 = float(np.percentile(cell_bacteria_raw, 98))
            signal_range = p98 - p2

            if signal_range < min_range:
                skipped_low_signal += 1
                empty_mask = np.zeros(bacteria_2d.shape, dtype=bool)
                results[cell_id] = {
                    'bacterial_area': 0,
                    'host_area': host_area,
                    'area_ratio': 0.0,
                    'bacteria_mask': empty_mask
                }
                continue

        # Extract adjusted ROI for detection
        bacteria_roi = bacteria_2d[r0:r1, c0:c1].copy()

        if bacteria_roi.size == 0:
            empty_mask = np.zeros(bacteria_2d.shape, dtype=bool)
            results[cell_id] = {
                'bacterial_area': 0,
                'host_area': host_area,
                'area_ratio': 0.0,
                'bacteria_mask': empty_mask
            }
            continue

        # Per-cell normalization
        if norm_params.get('enabled', True):
            bacteria_roi = normalize_cell_region(
                bacteria_roi,
                method=norm_params.get('method', 'percentile'),
                percentile_low=norm_params.get('percentile_low', 2),
                percentile_high=norm_params.get('percentile_high', 98)
            )

        # Sato detection
        bacteria_mask_roi = detect_bacteria_in_cell(bacteria_roi, detection_params)

        # Restrict to cell interior
        bacteria_mask_roi = bacteria_mask_roi & cell_mask_roi

        # Compute area
        bacterial_area = int(np.sum(bacteria_mask_roi))
        area_ratio = bacterial_area / float(host_area) if host_area > 0 else 0.0

        # Store full-image-size bacteria mask for zone analysis
        bacteria_mask_full = np.zeros(bacteria_2d.shape, dtype=bool)
        bacteria_mask_full[r0:r1, c0:c1] = bacteria_mask_roi

        results[cell_id] = {
            'bacterial_area': bacterial_area,
            'host_area': host_area,
            'area_ratio': area_ratio,
            'bacteria_mask': bacteria_mask_full
        }

    if skipped_low_signal > 0:
        print(f"  Skipped {skipped_low_signal}/{len(valid_cell_ids)} cells (low signal)")
    detected = sum(1 for r in results.values() if r['bacterial_area'] > 0)
    print(f"  Quantified {len(results)} cells ({detected} with bacteria)")

    return results


def compute_zone_bacterial_density(bacteria_mask, cell_mask, nucleus_mask, zone_thresholds=None):
    """
    Compute bacterial area and density within 3 nucleus-anchored zones:
    - Perinuclear (inner): normalized distance 0.0 - 0.33 from nucleus
    - Middle:              normalized distance 0.33 - 0.66
    - Cortical (outer):    normalized distance 0.66 - 1.0 (near membrane)

    Normalized distance = D_from_nucleus / (D_from_nucleus + D_from_membrane).
    At nucleus edge = 0; at cell membrane = 1.

    Parameters
    ----------
    bacteria_mask : ndarray [Y, X] bool
        Full-image-size binary bacteria detection mask for this cell.
        (Returned as 'bacteria_mask' key from quantify_bacteria_per_cell)
    cell_mask : ndarray [Y, X] bool
        Binary mask for this single cell.
    nucleus_mask : ndarray [Y, X] bool
        Binary mask for the nucleus belonging to this cell.
    zone_thresholds : dict or None
        Keys 'inner' and 'mid'. Defaults to {'inner': 0.33, 'mid': 0.66}.

    Returns
    -------
    dict with keys:
        bact_area_perinuclear, bact_area_middle, bact_area_cortical  (int, pixels)
        bact_density_perinuclear, bact_density_middle, bact_density_cortical
            (float, bacterial_pixels / zone_pixels)
    """
    from scipy import ndimage as _ndi

    if zone_thresholds is None:
        zone_thresholds = {'inner': 0.33, 'mid': 0.66}

    nan_result = {
        'bact_area_perinuclear':    0,
        'bact_area_middle':         0,
        'bact_area_cortical':       0,
        'bact_density_perinuclear': 0.0,
        'bact_density_middle':      0.0,
        'bact_density_cortical':    0.0,
    }

    if cell_mask is None or cell_mask.sum() == 0:
        return nan_result

    cell_int = cell_mask.astype(np.uint8)

    # Build normalized distance map
    if nucleus_mask is not None and nucleus_mask.sum() > 0:
        nuc_int = nucleus_mask.astype(np.uint8)
        dist_from_nuc = _ndi.distance_transform_edt(1 - nuc_int)
        dist_from_mem = _ndi.distance_transform_edt(cell_int)
        total_dist = dist_from_nuc + dist_from_mem
        with np.errstate(divide='ignore', invalid='ignore'):
            zone_map = np.divide(dist_from_nuc, total_dist,
                                 out=np.zeros_like(total_dist, dtype=np.float32),
                                 where=total_dist > 0)
    else:
        # Fallback: centroid-based (no nucleus available)
        dist_from_mem = _ndi.distance_transform_edt(cell_int)
        max_d = dist_from_mem.max()
        if max_d == 0:
            return nan_result
        zone_map = 1.0 - (dist_from_mem / max_d)

    # Zone masks (restricted to cell interior)
    t_inner = zone_thresholds['inner']
    t_mid   = zone_thresholds['mid']
    inner_mask = cell_mask & (zone_map >= 0.0) & (zone_map <  t_inner)
    mid_mask   = cell_mask & (zone_map >= t_inner) & (zone_map < t_mid)
    outer_mask = cell_mask & (zone_map >= t_mid)

    # Bacterial pixels per zone
    bact_inner = int(np.sum(bacteria_mask & inner_mask))
    bact_mid   = int(np.sum(bacteria_mask & mid_mask))
    bact_outer = int(np.sum(bacteria_mask & outer_mask))

    # Zone areas
    area_inner = int(inner_mask.sum())
    area_mid   = int(mid_mask.sum())
    area_outer = int(outer_mask.sum())

    return {
        'bact_area_perinuclear':    bact_inner,
        'bact_area_middle':         bact_mid,
        'bact_area_cortical':       bact_outer,
        'bact_density_perinuclear': round(bact_inner / area_inner, 6) if area_inner > 0 else 0.0,
        'bact_density_middle':      round(bact_mid   / area_mid,   6) if area_mid   > 0 else 0.0,
        'bact_density_cortical':    round(bact_outer / area_outer, 6) if area_outer > 0 else 0.0,
    }
