"""Per-cell Sato-weighted skeleton density, as used for the figure.

  1. Gaussian sigma 0.7, then the multiscale Sato filter (sigmas 0.7-2.0, 6 steps,
     scikit-image defaults) on the whole field. The input is the actin maximum projection
     rescaled to 0-1 (local_adjust_image).
  2. Sato rescaled 0-1 over the cell's pixels; threshold = 55th percentile of the positive
     values; objects < 3 px removed; closing with a disk of radius 1; skeletonize.
  3. Radial position r = d_nucleus / (d_nucleus + d_membrane) relative to the cell's own
     nucleus (0 at the nucleus, 1 at the membrane).
  4. Sato-weighted skeleton density = (Sato response summed on the skeleton) / (pixel area),
     per radial bin (10 bins of r) and per zone (inner r < 0.33, mid 0.33-0.66, outer >= 0.66),
     and for the ten rings (ring 0 = nucleus, rings 1-9 = cytoplasm).
"""
import numpy as np
from scipy import ndimage
from skimage import filters, morphology
from skimage.morphology import disk, skeletonize

from v3_paper import ring_labels, ring_profile

SATO_SIGMAS = [0.7, 0.96, 1.22, 1.48, 1.74, 2.0]


def local_adjust_image(img, bright=0.0, contrast=1.0, gamma=1.0):
    """Rescale to 0-1, then apply gamma, contrast and brightness (defaults: no change)."""
    img = img.astype(float)
    x = (img - img.min()) / (img.max() - img.min() + 1e-9)
    x = np.power(x, gamma)
    x = (x - 0.5) * contrast + 0.5
    return np.clip(x + bright, 0, 1)


def sato_field(actin_2d, sigmas=SATO_SIGMAS, black_ridges=True):
    smoothed = ndimage.gaussian_filter(actin_2d.astype(float), sigma=0.7)
    return filters.sato(smoothed, sigmas=sigmas, mode="reflect", black_ridges=black_ridges)


def skeleton_of_cell(sato_full, cell_mask_full, percentile=55.0, min_size=3, closing=1):
    inside = sato_full[cell_mask_full]
    lo, hi = inside.min(), inside.max()
    norm = (sato_full - lo) / (hi - lo) if hi > lo else sato_full.copy()
    cell = norm * cell_mask_full
    pos = cell[cell_mask_full]
    pos = pos[pos > 0]
    thr = np.percentile(pos, percentile) if len(pos) > 10 else 0.5
    binary = (cell > thr) & cell_mask_full
    binary = morphology.remove_small_objects(binary, min_size=min_size)
    if closing:
        binary = morphology.binary_closing(binary, disk(closing))
    return binary, skeletonize(binary)


def normalized_distance(cell_mask, nucleus_mask):
    """V3's TopologicalAnalyzer.get_normalized_distance_map (nucleus-anchored case)."""
    d_nuc = ndimage.distance_transform_edt(1 - nucleus_mask)
    d_mem = ndimage.distance_transform_edt(cell_mask)
    total = d_nuc + d_mem
    r = np.divide(d_nuc, total, out=np.zeros_like(total), where=total != 0)
    return np.where(cell_mask, r, np.nan)


def bin_profile(skeleton, sato, r, cell_mask, n_bins=10):
    """Density in n_bins equal steps of r from 0 to 1 (the nucleus, r = 0, falls in the first bin)."""
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    safe = np.nan_to_num(r, nan=-1.0)
    dens = np.zeros(n_bins)
    for i in range(n_bins):
        upper = (safe <= edges[i + 1]) if i == n_bins - 1 else (safe < edges[i + 1])
        m = (safe >= edges[i]) & upper & cell_mask
        if m.sum():
            dens[i] = sato[skeleton & m].sum() / m.sum()
    return (edges[:-1] + edges[1:]) / 2, dens


def zone_densities(skeleton, sato, r, cell_mask):
    safe = np.nan_to_num(r, nan=-1.0)
    zones = {"inner": (safe >= 0.0) & (safe < 0.33), "mid": (safe >= 0.33) & (safe < 0.66),
             "cortical": safe >= 0.66}
    out = {}
    for name, z in zones.items():
        z = z & cell_mask
        out[f"{name}_weighted_density"] = sato[skeleton & z].sum() / z.sum() if z.sum() else np.nan
    return out


def analyze_cell(actin_2d, cell_masks, nucleus_masks, cell_id, nucleus_id, sato_full=None,
                 black_ridges=True, percentile=55.0, n_bins=10):
    """Measure one cell. Returns the crops, skeleton, Sato response, r map and the densities.

    actin_2d is the image the Sato filter runs on: local_adjust_image(actin maximum projection).
    """
    if sato_full is None:
        sato_full = sato_field(actin_2d, black_ridges=black_ridges)
    cell_full = cell_masks == cell_id
    bbox = ndimage.find_objects(cell_full.astype(int))[0]
    mask = cell_full[bbox]
    nucleus = (nucleus_masks == nucleus_id)[bbox] if nucleus_id is not None else np.zeros_like(mask)
    binary_full, skeleton_full = skeleton_of_cell(sato_full, cell_full, percentile=percentile)
    binary, skeleton = binary_full[bbox], skeleton_full[bbox] & mask
    sato = sato_full[bbox] * mask
    r = normalized_distance(mask, nucleus)
    centres, bins = bin_profile(skeleton, sato, r, mask, n_bins=n_bins)
    row = dict(zone_densities(skeleton, sato, r, mask))
    row.update(ring_profile(skeleton, sato, ring_labels(r, mask)))
    return {"cell_id": cell_id, "bbox": bbox, "mask": mask, "nucleus": nucleus, "binary": binary, "skeleton": skeleton,
            "sato": sato, "norm_map": r, "rings": ring_labels(r, mask),
            "bin_centres": centres, "bin_profile": bins, "row": row}
