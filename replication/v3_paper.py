"""Paper settings for the old TopoActin V3 code, and single-field analysis.

The code in old_pipeline/ is kept exactly as it was. This module applies the
two run-time switches the paper runs need, adds the 10-ring profile, and runs
one field step by step with the same V3 classes that main_pipeline_v3.py uses,
so that a notebook can show every intermediate image.

    from v3_paper import load_config, apply_paper_settings, analyze_field
    cfg = load_config()
    apply_paper_settings(cfg)
    field = analyze_field("C9_F13_Z_stack_merged.tiff", cfg)
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
OLD = REPO / "old_pipeline"
if str(OLD) not in sys.path:
    sys.path.insert(0, str(OLD))

import segmentation_v3  # noqa: E402
import skeletonization_v3  # noqa: E402
import topological_analysis_v3  # noqa: E402
from filtering import CellFilter, NucleiFilter  # noqa: E402

DEFAULT_CONFIG = REPO / "replication" / "config_paper.json"
N_RINGS = 9  # cytoplasm rings; ring 0 is the nucleus
ZONES = {"inner": (1, 2, 3), "mid": (4, 5, 6), "outer": (7, 8, 9)}

_applied = set()


def load_config(path=DEFAULT_CONFIG, cpu=False):
    cfg = json.loads(Path(path).read_text(encoding="utf-8"))
    if cpu:
        cfg.setdefault("processing", {})["use_gpu"] = False
    return cfg


def apply_paper_settings(cfg, measure_on="adjusted"):
    """Apply the two code-level switches the January 2026 paper runs had.

    measure_on="adjusted": intensities are measured on the brightness/contrast-
    adjusted actin image. The V3 code was changed on 10 Feb 2026 to measure on the
    raw image, which gives values 12-25x lower than the paper.

    sato_black_ridges=false in the config: skimage's sato() defaults to dark
    ridges, which traces the gaps between fibres; actin fibres are bright.
    """
    if measure_on == "adjusted" and "adjusted" not in _applied:
        original = segmentation_v3.CellSegmenter.segment_image

        def segment_and_measure_adjusted(self, path):
            result = original(self, path)
            if result is not None:
                result["actin_image_unadjusted"] = result.get("actin_image_raw")
                result["actin_image_raw"] = result["actin_image"]
                if "actin_volume_3d" in result:
                    result["actin_volume_3d_raw"] = result["actin_volume_3d"]
            return result

        segmentation_v3.CellSegmenter.segment_image = segment_and_measure_adjusted
        _applied.add("adjusted")

    if not cfg.get("skeletonization", {}).get("sato_black_ridges", True) and "bright" not in _applied:
        from skimage import filters as _skf

        def sato_bright(self, volume):
            return _skf.sato(volume, sigmas=self.sato_sigmas, mode="reflect", black_ridges=False)

        skeletonization_v3.CellSkeletonizer._apply_sato_3d = sato_bright
        _applied.add("bright")


def ring_labels(norm_map, cell_mask):
    """Label image of the 10 rings: 0 = nucleus, 1..9 = cytoplasm rings, -1 = outside.

    norm_map is V3's nucleus-anchored distance r = d_nucleus / (d_nucleus + d_membrane):
    0 inside the nucleus and at its edge, 1 at the membrane. Ring k (1..9) holds
    (k-1)/9 < r <= k/9.
    """
    cell = cell_mask.astype(bool)
    r = np.where(cell, norm_map, np.nan)
    labels = np.full(cell.shape, -1, dtype=np.int8)
    labels[cell & (r == 0)] = 0
    for k in range(1, N_RINGS + 1):
        labels[cell & (r > (k - 1) / N_RINGS) & (r <= k / N_RINGS)] = k
    return labels


def ring_profile(skeleton, sato_map, labels):
    """Sato-weighted skeleton density per ring and per merged 3-3-3 zone.

    density = (Sato summed on the skeleton pixels in the ring) / (ring area in pixels),
    the formula of V3's calculate_weighted_skeleton_density. Merged zones sum the
    signal and the area of their rings (they are not averages of ring densities).
    """
    skel = skeleton.astype(bool)
    signal = np.array([float(np.sum(sato_map[skel & (labels == k)])) for k in range(N_RINGS + 1)])
    area = np.array([int(np.sum(labels == k)) for k in range(N_RINGS + 1)])
    out = {}
    with np.errstate(invalid="ignore", divide="ignore"):
        for k in range(N_RINGS + 1):
            out[f"ring{k}_weighted_density"] = signal[k] / area[k] if area[k] else np.nan
        for name, ks in ZONES.items():
            a = area[list(ks)].sum()
            out[f"ring_{name}3_weighted_density"] = signal[list(ks)].sum() / a if a else np.nan
    return out


def add_ring_profile():
    """Make the V3 batch run also write the 10-ring profile columns for every cell.

    Ring 0 is the nucleus, rings 1-9 the cytoplasm; ring_inner3 / ring_mid3 /
    ring_outer3 are rings 1-3 / 4-6 / 7-9 merged. V3's own zones use the cut-offs
    0.33 / 0.66 and include the nucleus pixels (r = 0) in the inner zone, so
    ring_outer3 ~ cortical_weighted_density, while inner_weighted_density also
    contains the nucleus.
    """
    if "rings" in _applied:
        return
    last = {}
    get_map = topological_analysis_v3.TopologicalAnalyzer.get_normalized_distance_map

    def get_map_and_keep(self, cell_mask, nucleus_mask):
        norm_map, method = get_map(self, cell_mask, nucleus_mask)
        last["norm"], last["cell"], last["method"] = norm_map, cell_mask, method
        return norm_map, method

    topological_analysis_v3.TopologicalAnalyzer.get_normalized_distance_map = get_map_and_keep
    ratios = skeletonization_v3.calculate_weighted_skeleton_ratios

    def ratios_and_rings(skeleton, sato_map, zones):
        out = ratios(skeleton, sato_map, zones)
        norm = last.get("norm")
        if norm is not None and last.get("method") == "nucleus_anchored" and norm.shape == skeleton.shape:
            out.update(ring_profile(skeleton & last["cell"].astype(bool), sato_map,
                                    ring_labels(norm, last["cell"])))
        return out

    skeletonization_v3.calculate_weighted_skeleton_ratios = ratios_and_rings
    _applied.add("rings")


def analyze_field(image_path, cfg, verbose=False, skeletons=True):
    """Run V3 on one field and keep every intermediate result.

    Same order of calls as main_pipeline_v3.process_experiment: Cellpose
    segmentation, nucleus filter, cell filter, then skeletonize_cell per kept cell.
    Returns a dict with the images, masks, kept cells, per-cell V3 results
    (skeleton, Sato map, bounding box, features) and the 10-ring profile.
    skeletons=False stops after the filters (segmentation and kept cells only).
    """
    import contextlib
    import io

    cfg = dict(cfg)
    cfg.setdefault("output_dir", str(Path(image_path).resolve().parent))
    cfg.setdefault("processing", {})["save_skeletons"] = False
    cfg["processing"]["save_sato_maps"] = False
    quiet = contextlib.nullcontext() if verbose else contextlib.redirect_stdout(io.StringIO())

    with quiet:
        segmenter = segmentation_v3.CellSegmenter(cfg)
        seg = segmenter.segment_image(str(image_path))
        n_res = NucleiFilter(cfg).filter_nuclei(seg["nucleus_masks"], seg["nuclei_image"], Path(image_path).name)
        c_res = CellFilter(cfg).filter_cells(seg["cell_masks"], seg["nucleus_masks"],
                                             n_res["valid_nuclei"], Path(image_path).name)
        cells = []
        if not skeletons:
            return {"image_path": str(image_path), "seg": seg, "nuclei_filter": n_res,
                    "cell_filter": c_res, "cells": cells}
        skeletonizer = skeletonization_v3.CellSkeletonizer(cfg)
        analyzer = topological_analysis_v3.TopologicalAnalyzer(cfg.get("topological_analysis", {}))
        for info in c_res["valid_cells"]:
            cid = info["cell_id"]
            cell_mask = seg["cell_masks"] == cid
            res = skeletonizer.skeletonize_cell(
                cell_mask=cell_mask, actin_image=seg["actin_image"], image_name=Path(image_path).name,
                cell_id=cid, nucleus_mask=seg["nucleus_masks"],
                actin_volume_3d=seg.get("actin_volume_3d"),
                actin_image_raw=seg.get("actin_image_raw", seg["actin_image"]),
                actin_volume_3d_raw=seg.get("actin_volume_3d_raw", seg.get("actin_volume_3d")))
            if not res:
                continue
            bbox = res["bbox"]
            mask_crop = cell_mask[bbox]
            # as in skeletonize_cell: every nucleus pixel inside the bounding box
            nuc_crop = seg["nucleus_masks"][bbox] > 0
            norm, method = analyzer.get_normalized_distance_map(mask_crop, nuc_crop)
            labels = ring_labels(norm, mask_crop)
            row = {"cell_id": cid, **info["properties"], **res["features"], **res["cortical_features"]}
            if method == "nucleus_anchored":
                row.update(ring_profile(res["skeleton"] & mask_crop, res["sato_raw"], labels))
            cells.append({"cell_id": cid, "bbox": bbox, "mask": mask_crop, "nucleus": nuc_crop,
                          "skeleton": res["skeleton"], "sato": res["sato_raw"], "norm_map": norm,
                          "rings": labels, "row": row})
    return {"image_path": str(image_path), "seg": seg, "nuclei_filter": n_res, "cell_filter": c_res,
            "cells": cells}


def cell_table(field):
    import pandas as pd
    return pd.DataFrame([c["row"] for c in field["cells"]])
