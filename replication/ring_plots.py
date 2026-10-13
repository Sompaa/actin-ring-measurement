"""Figures for notebooks/single_image_walkthrough.ipynb.

Style of the figure: actin red, nuclei blue; zones blue (inner, r < 0.33),
yellow (mid, 0.33-0.66), red (outer / cortical, >= 0.66); skeleton white; cell outline cyan;
radial profile as a blue line with points over a light-blue fill, dashed lines at 0.33 / 0.66.
"""
import matplotlib.pyplot as plt
import numpy as np
from skimage.segmentation import find_boundaries

ZONE_RGB = {"inner": (0.2, 0.4, 0.9), "mid": (0.9, 0.8, 0.2), "outer": (0.9, 0.3, 0.2)}


def _unit(img):
    img = img.astype(np.float32)
    return (img - img.min()) / (img.max() - img.min() + 1e-9)


def composite(actin, nuclei):
    """Red actin + blue nuclei, each rescaled 0-1 over the field (as in the notebook)."""
    rgb = np.zeros(actin.shape + (3,), np.float32)
    rgb[..., 0] = _unit(actin)
    rgb[..., 2] = _unit(nuclei)
    return rgb


def zone_masks(r, mask):
    safe = np.nan_to_num(r, nan=-1.0)
    return {"inner": (safe >= 0) & (safe < 0.33) & mask, "mid": (safe >= 0.33) & (safe < 0.66) & mask,
            "outer": (safe >= 0.66) & mask}


def show_channels(actin, nuclei):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(composite(actin, nuclei))
    ax.set_title("Actin (red) and nuclei (blue), max projection")
    ax.axis("off")
    return fig


def show_segmentation(actin, nuclei, cell_masks, nucleus_masks, kept):
    """All Cellpose cells (grey outline), the cells kept by the filters (cyan, numbered)."""
    rgb = composite(actin, nuclei)
    rgb[find_boundaries(cell_masks, mode="inner")] = (0.6, 0.6, 0.6)
    kept_lab = np.where(np.isin(cell_masks, list(kept)), cell_masks, 0)
    rgb[find_boundaries(kept_lab, mode="inner")] = (0, 1, 1)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(rgb)
    for cid in kept:
        ys, xs = np.nonzero(cell_masks == cid)
        ax.text(xs.mean(), ys.mean(), str(cid), color="white", ha="center", va="center", fontsize=9,
                fontweight="bold")
    n_all = len(np.unique(cell_masks)) - 1
    ax.set_title(f"Cellpose found {n_all} cells (grey); {len(kept)} pass the filters (cyan)")
    ax.axis("off")
    return fig


def eight_panels(cell, image):
    """Intermediate steps for one cell, eight panels."""
    bbox, mask, skel = cell["bbox"], cell["mask"], cell["skeleton"]
    crop = np.where(mask, _unit(image[bbox]), 0)
    fig, ax = plt.subplots(2, 4, figsize=(20, 9))
    ax = ax.ravel()
    ax[0].imshow(crop, cmap="gray")
    ax[0].contour(mask, levels=[0.5], colors="cyan", linewidths=1)
    ax[0].set_title("1. Original cell region")
    im = ax[1].imshow(cell["sato"], cmap="hot")
    fig.colorbar(im, ax=ax[1], fraction=0.046)
    ax[1].set_title("2. Sato tubeness")
    ax[2].imshow(cell["binary"] & mask, cmap="gray")
    ax[2].set_title("3. Binary threshold (55th percentile)")
    g = np.dstack([crop * 0.6 + 0.15] * 3)
    g[skel] = (0, 1, 1)
    ax[3].imshow(g)
    ax[3].set_title("4. Final skeleton")
    t = np.dstack([_unit(cell["sato"])] * 3) * np.array([1.0, 0.25, 0.1])
    t[skel] = 1
    ax[4].imshow(t)
    ax[4].set_title("5. Tubeness + skeleton")
    z = np.zeros(mask.shape + (3,))
    for name, m in zone_masks(cell["norm_map"], mask).items():
        z[m] = ZONE_RGB[name]
    z[skel] = 1
    ax[5].imshow(z)
    ax[5].set_title("6. Radial zones (blue = inner, yellow = mid, red = outer)")
    radial_profile(cell, ax=ax[6])
    ax[6].set_title("7. Radial profile")
    row = cell["row"]
    ax[7].axis("off")
    ax[7].text(0.02, 0.95, (
        f"Cell {cell['cell_id']}\n\nSato-weighted skeleton density\n"
        f"  inner (0-0.33):     {row['inner_weighted_density']:.4f}\n"
        f"  mid   (0.33-0.66):  {row['mid_weighted_density']:.4f}\n"
        f"  outer (0.66-1):     {row['cortical_weighted_density']:.4f}\n\n"
        f"Outer / inner: {row['cortical_weighted_density'] / row['inner_weighted_density']:.2f}"),
        va="top", family="monospace", fontsize=11,
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.6))
    ax[7].set_title("8. Metrics")
    for a in ax[[0, 1, 2, 3, 4, 5]]:
        a.axis("off")
    fig.suptitle(f"Cytoskeleton analysis - cell {cell['cell_id']}", fontweight="bold")
    fig.tight_layout()
    return fig


def radial_profile(cell, ax=None):
    """Sato-weighted skeleton density in bins of r (0 = nucleus, 1 = membrane)."""
    fig, ax = (plt.subplots(figsize=(5, 4)) if ax is None else (ax.figure, ax))
    x, y = cell["bin_centres"], cell["bin_profile"]
    ax.fill_between(x, 0, y, color="steelblue", alpha=0.3)
    ax.plot(x, y, "o-", color="blue", label="Weighted skeleton density")
    ax.axvline(0.33, color="orange", ls="--", label="Inner/Mid")
    ax.axvline(0.66, color="red", ls="--", label="Mid/Outer")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, None)
    ax.set_xlabel("Radial distance (0 = nucleus, 1 = membrane)")
    ax.set_ylabel("Weighted skeleton density")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    return fig


def paper_figure(cell, actin, nuclei, cell_masks, alpha=0.45, margin=150):
    """The three panels of the figure: original image, radial overlay of the cell
    (zones, skeleton, outline) and the radial distribution plot.

    The two image panels show the cell's bounding box plus `margin` px on each side; margin=None shows the whole field.
    """
    rgb = composite(actin, nuclei)
    over = rgb.copy()
    bbox, mask = cell["bbox"], cell["mask"]
    sub = over[bbox]
    for name, m in zone_masks(cell["norm_map"], mask).items():
        sub[m] = sub[m] * (1 - alpha) + np.array(ZONE_RGB[name]) * alpha
    sub[cell["skeleton"]] = 1
    full_mask = cell_masks == cell["cell_id"]
    if margin is not None:
        win = tuple(slice(max(s.start - margin, 0), min(s.stop + margin, n))
                    for s, n in zip(bbox, full_mask.shape))
        rgb, over, full_mask = rgb[win], over[win], full_mask[win]
    fig, ax = plt.subplots(1, 3, figsize=(19, 5.5), gridspec_kw={"width_ratios": [1, 1, 0.9]})
    ax[0].imshow(rgb)
    ax[0].set_title("Original image")
    ax[1].imshow(over)
    ax[1].contour(full_mask, levels=[0.5], colors="cyan", linewidths=1)
    ax[1].set_title("Radial overlay of cell")
    for a in ax[:2]:
        a.axis("off")
    radial_profile(cell, ax=ax[2])
    ax[2].set_title("Distribution plot")
    fig.tight_layout()
    return fig


def ten_rings(cell, ax=None):
    """Ring 0 = nucleus, rings 1-9 = cytoplasm, and the rings merged 3-3-3."""
    row = cell["row"]
    fig, ax = (plt.subplots(figsize=(7, 4)) if ax is None else (ax.figure, ax))
    x = np.arange(10)
    y = [row[f"ring{k}_weighted_density"] for k in x]
    colours = ["0.6"] + [ZONE_RGB["inner"]] * 3 + [ZONE_RGB["mid"]] * 3 + [ZONE_RGB["outer"]] * 3
    ax.bar(x, y, color=colours, edgecolor="k", lw=0.5)
    for name, ks in (("inner", (1, 3)), ("mid", (4, 6)), ("outer", (7, 9))):
        v = row[f"ring_{name}3_weighted_density"]
        ax.hlines(v, ks[0] - 0.45, ks[1] + 0.45, color="k", lw=2.5)
        ax.text((ks[0] + ks[1]) / 2, v, f"{v:.4f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x, ["nucleus"] + [str(k) for k in range(1, 10)])
    ax.set_xlabel("ring")
    ax.set_ylabel("Sato-weighted skeleton density")
    ax.set_title("10 rings (black bars: rings merged 3-3-3)")
    return fig
