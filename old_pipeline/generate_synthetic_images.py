"""
generate_synthetic_images.py — Create synthetic 4-channel TIFF test images for
pipeline accuracy validation.

Usage:
    python generate_synthetic_images.py [output_dir]

Default output_dir: synthetic_test/ next to this script.

Produces:
    synthetic_test/
        stress_fibers.tif       -- parallel linear actin fibers
        cortical_ring.tif       -- peripheral actin ring, sparse center
        dense_mesh.tif          -- random actin filament network
        infected_cell.tif       -- moderate actin + bacterial rods near nucleus
        empty_control.tif       -- near-background actin signal
        ground_truth.csv        -- known metrics for accuracy comparison

Image format: [Z=5, C=4, Y=512, X=512] uint16
Channel order: actin=0, nuclei=1, brightfield=2, bacteria=3
"""

import sys
import csv
import pathlib
import numpy as np
from scipy.ndimage import gaussian_filter

try:
    import tifffile
except ImportError:
    print("ERROR: tifffile not installed. Run: pip install tifffile")
    sys.exit(1)

try:
    from skimage import draw, morphology, filters
    from skimage.morphology import disk, dilation, binary_dilation
except ImportError:
    print("ERROR: scikit-image not installed. Run: pip install scikit-image")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
H, W = 512, 512
N_Z = 5
N_C = 4                   # actin, nuclei, brightfield, bacteria
DTYPE = np.uint16
MAX_VAL = 65535

# Intensity levels (uint16)
BG_MEAN   = 500
BG_SIGMA  = 200
ACTIN_SIG = 18000         # filament peak signal
NUCLEUS_SIG = 10000
BRIGHTFIELD_SIG = 30000
BACTERIA_SIG = 16000

# Cell geometry (centred in image)
CY, CX = H // 2, W // 2
CELL_RY, CELL_RX = 180, 140   # cell semi-axes (px)
NUC_R = 50                     # nucleus radius (px)


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def make_cell_mask() -> np.ndarray:
    """Boolean mask of elliptical cell body."""
    rr, cc = draw.ellipse(CY, CX, CELL_RY, CELL_RX, shape=(H, W))
    mask = np.zeros((H, W), dtype=bool)
    mask[rr, cc] = True
    return mask


def make_nucleus_mask() -> np.ndarray:
    """Boolean mask of circular nucleus."""
    rr, cc = draw.disk((CY, CX), NUC_R, shape=(H, W))
    mask = np.zeros((H, W), dtype=bool)
    mask[rr, cc] = True
    return mask


def make_brightfield(cell_mask: np.ndarray) -> np.ndarray:
    """Simple brightfield: bright background, dark cell interior."""
    bf = np.full((H, W), BRIGHTFIELD_SIG, dtype=np.float32)
    bf[cell_mask] = BRIGHTFIELD_SIG * 0.55
    bf = gaussian_filter(bf, sigma=3.0)
    return np.clip(bf, 0, MAX_VAL).astype(DTYPE)


def make_nucleus_channel(nuc_mask: np.ndarray) -> np.ndarray:
    """Nuclear staining channel."""
    ch = np.random.normal(BG_MEAN, BG_SIGMA, (H, W)).astype(np.float32)
    nuc_float = nuc_mask.astype(np.float32)
    nuc_smooth = gaussian_filter(nuc_float, sigma=4.0)
    ch += nuc_smooth * NUCLEUS_SIG
    return np.clip(ch, 0, MAX_VAL).astype(DTYPE)


def add_noise(img: np.ndarray) -> np.ndarray:
    noise = np.random.normal(0, BG_SIGMA, img.shape).astype(np.float32)
    return np.clip(img.astype(np.float32) + noise, 0, MAX_VAL).astype(DTYPE)


def apply_defocus(img: np.ndarray, z: int) -> np.ndarray:
    """Simulate defocus for non-central Z slices."""
    sigma = abs(z - N_Z // 2) * 0.8
    if sigma > 0:
        img = gaussian_filter(img.astype(np.float32), sigma=sigma)
    return img.astype(DTYPE)


def build_stack(actin_2d: np.ndarray, nuc_2d: np.ndarray,
                bf_2d: np.ndarray, bact_2d: np.ndarray) -> np.ndarray:
    """
    Build [Z, C, Y, X] stack with defocus variation across Z.
    Middle slice (z=2) is sharpest so MIP picks it up correctly.
    """
    stack = np.zeros((N_Z, N_C, H, W), dtype=DTYPE)
    for z in range(N_Z):
        stack[z, 0] = add_noise(apply_defocus(actin_2d, z))
        stack[z, 1] = add_noise(apply_defocus(nuc_2d, z))
        stack[z, 2] = apply_defocus(bf_2d, z)
        stack[z, 3] = add_noise(apply_defocus(bact_2d, z))
    return stack


def empty_bacteria() -> np.ndarray:
    return np.zeros((H, W), dtype=DTYPE)


# ---------------------------------------------------------------------------
# Actin pattern generators
# ---------------------------------------------------------------------------

def draw_fiber(canvas: np.ndarray, r0, c0, r1, c1, width: int = 3) -> int:
    """Draw a single fiber and return its pixel length."""
    rr, cc = draw.line(r0, c0, r1, c1)
    skel = np.zeros((H, W), dtype=bool)
    skel[rr, cc] = True
    thick = dilation(skel, disk(width))
    canvas[thick] = ACTIN_SIG
    return len(rr)


def make_stress_fibers(cell_mask: np.ndarray) -> tuple:
    """
    ~10 parallel fibers running horizontally through the cell.
    Returns (actin_channel, ground_truth_dict).
    """
    actin = np.zeros((H, W), dtype=np.float32)
    n_fibers = 10
    fiber_spacing = (CELL_RY * 2) // (n_fibers + 1)
    total_length = 0
    fiber_width = 2

    for i in range(1, n_fibers + 1):
        fy = CY - CELL_RY + i * fiber_spacing
        # clip to within cell horizontally
        fc0 = CX - CELL_RX + 15
        fc1 = CX + CELL_RX - 15
        length = draw_fiber(actin, fy, fc0, fy, fc1, width=fiber_width)
        total_length += length

    actin = gaussian_filter(actin, sigma=1.2)
    actin[~cell_mask] = 0

    gt = {
        "scenario": "stress_fibers",
        "approx_total_skeleton_length_px": total_length,
        "approx_num_fibers": n_fibers,
        "expected_tortuosity": "~1.0 (straight fibers)",
        "expected_cortical_vs_inner": "lower cortical (fibers span full cell)",
        "has_bacteria": False,
    }
    return np.clip(actin, 0, MAX_VAL).astype(DTYPE), gt


def make_cortical_ring(cell_mask: np.ndarray) -> tuple:
    """
    Thick actin ring at cell periphery, sparse interior.
    Returns (actin_channel, ground_truth_dict).
    """
    actin = np.zeros((H, W), dtype=np.float32)

    # Draw ellipse outline (cortical ring)
    ring_skel = np.zeros((H, W), dtype=bool)
    rr, cc = draw.ellipse_perimeter(CY, CX, CELL_RY - 10, CELL_RX - 10, shape=(H, W))
    ring_skel[rr, cc] = True
    ring_thick = dilation(ring_skel, disk(12))
    actin[ring_thick & cell_mask] = ACTIN_SIG

    # Sparse filaments in interior
    for _ in range(15):
        angle = np.random.uniform(0, 2 * np.pi)
        r_start = int(CY + 20 * np.sin(angle))
        c_start = int(CX + 20 * np.cos(angle))
        r_end = int(CY + (NUC_R + 40) * np.sin(angle))
        c_end = int(CX + (NUC_R + 40) * np.cos(angle))
        draw_fiber(actin, r_start, c_start, r_end, c_end, width=1)

    actin = gaussian_filter(actin, sigma=1.5)
    actin[~cell_mask] = 0

    gt = {
        "scenario": "cortical_ring",
        "expected_cortical_density": "HIGH",
        "expected_inner_density": "LOW",
        "expected_cortical_inner_ratio": "> 2.0",
        "has_bacteria": False,
    }
    return np.clip(actin, 0, MAX_VAL).astype(DTYPE), gt


def make_dense_mesh(cell_mask: np.ndarray) -> tuple:
    """
    Random 2D actin filament network (lamellipodia-like).
    Returns (actin_channel, ground_truth_dict).
    """
    actin = np.zeros((H, W), dtype=np.float32)
    n_fibers = 60
    branch_points_approx = 0

    for _ in range(n_fibers):
        angle = np.random.uniform(0, 2 * np.pi)
        length = np.random.randint(40, 120)
        # Random start within cell
        r0 = np.random.randint(CY - CELL_RY + 20, CY + CELL_RY - 20)
        c0 = np.random.randint(CX - CELL_RX + 20, CX + CELL_RX - 20)
        r1 = int(r0 + length * np.sin(angle))
        c1 = int(c0 + length * np.cos(angle))
        r1 = np.clip(r1, 0, H - 1)
        c1 = np.clip(c1, 0, W - 1)
        draw_fiber(actin, r0, c0, r1, c1, width=1)

    # Rough estimate of branch points: more fibers → more crossings
    branch_points_approx = n_fibers * 3  # empirical rough estimate

    actin = gaussian_filter(actin, sigma=1.0)
    actin[~cell_mask] = 0

    gt = {
        "scenario": "dense_mesh",
        "approx_num_source_fibers": n_fibers,
        "approx_branch_points": branch_points_approx,
        "expected_num_branches": "HIGH (>50)",
        "expected_tortuosity": "1.1 - 1.5",
        "has_bacteria": False,
    }
    return np.clip(actin, 0, MAX_VAL).astype(DTYPE), gt


def make_infected_cell(cell_mask: np.ndarray, nuc_mask: np.ndarray) -> tuple:
    """
    Moderate actin network + 4 bacterial rods positioned near nucleus.
    Returns (actin_channel, bacteria_channel, ground_truth_dict).
    """
    # Moderate actin (mix of ring + a few fibers)
    actin = np.zeros((H, W), dtype=np.float32)

    # Cortical ring (thinner)
    ring_skel = np.zeros((H, W), dtype=bool)
    rr, cc = draw.ellipse_perimeter(CY, CX, CELL_RY - 10, CELL_RX - 10, shape=(H, W))
    ring_skel[rr, cc] = True
    ring_thick = dilation(ring_skel, disk(7))
    actin[ring_thick & cell_mask] = ACTIN_SIG * 0.7

    # A few radial fibers
    for i in range(6):
        angle = i * (2 * np.pi / 6)
        r0 = int(CY + (NUC_R + 5) * np.sin(angle))
        c0 = int(CX + (NUC_R + 5) * np.cos(angle))
        r1 = int(CY + (CELL_RY - 20) * np.sin(angle))
        c1 = int(CX + (CELL_RX - 20) * np.cos(angle))
        draw_fiber(actin, r0, c0, r1, c1, width=2)

    actin = gaussian_filter(actin, sigma=1.3)
    actin[~cell_mask] = 0

    # Bacteria: 4 rods near nucleus (perinuclear zone = r < NUC_R * 2.5)
    bacteria = np.zeros((H, W), dtype=np.float32)
    bact_positions = []
    rod_length = 18
    for i in range(4):
        angle = i * (2 * np.pi / 4) + 0.3
        br = int(CY + (NUC_R + 20) * np.sin(angle))
        bc = int(CX + (NUC_R + 20) * np.cos(angle))
        # Rod oriented along angle
        r0 = int(br - (rod_length // 2) * np.sin(angle))
        c0 = int(bc - (rod_length // 2) * np.cos(angle))
        r1 = int(br + (rod_length // 2) * np.sin(angle))
        c1 = int(bc + (rod_length // 2) * np.cos(angle))
        r0, r1 = np.clip([r0, r1], 0, H - 1)
        c0, c1 = np.clip([c0, c1], 0, W - 1)
        rod_skel = np.zeros((H, W), dtype=bool)
        rr, cc = draw.line(r0, c0, r1, c1)
        rod_skel[rr, cc] = True
        rod_thick = dilation(rod_skel, disk(3))
        bacteria[rod_thick] = BACTERIA_SIG
        bact_positions.append((br, bc))

    bacteria = gaussian_filter(bacteria, sigma=1.0)
    bacteria[~cell_mask] = 0

    bact_area = int(np.sum(bacteria > BACTERIA_SIG * 0.3))

    gt = {
        "scenario": "infected_cell",
        "approx_bacteria_count": 4,
        "approx_bact_pixel_area": bact_area,
        "expected_bact_density_perinuclear": "> 0",
        "expected_bact_density_cortical": "~0 (bacteria near nucleus)",
        "has_bacteria": True,
    }
    return (np.clip(actin, 0, MAX_VAL).astype(DTYPE),
            np.clip(bacteria, 0, MAX_VAL).astype(DTYPE),
            gt)


def make_empty_control(cell_mask: np.ndarray) -> tuple:
    """
    Near-background actin signal. Tests pipeline robustness on low-signal cells.
    """
    actin = np.zeros((H, W), dtype=np.float32)
    # Very faint diffuse signal (unhealthy / low-expression cell)
    actin[cell_mask] = BG_MEAN * 0.8
    actin = gaussian_filter(actin, sigma=3.0)

    gt = {
        "scenario": "empty_control",
        "expected_skeleton_density": "~0 or very low",
        "expected_pipeline_behaviour": "Should not crash; may return 0-metric row",
        "has_bacteria": False,
    }
    return np.clip(actin, 0, MAX_VAL).astype(DTYPE), gt


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(out_dir: pathlib.Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    np.random.seed(42)

    cell_mask = make_cell_mask()
    nuc_mask  = make_nucleus_mask()
    nuc_ch    = make_nucleus_channel(nuc_mask)
    bf_ch     = make_brightfield(cell_mask)
    empty_bact = empty_bacteria()

    ground_truth_rows = []

    scenarios = [
        ("stress_fibers",  None),
        ("cortical_ring",  None),
        ("dense_mesh",     None),
        ("infected_cell",  None),
        ("empty_control",  None),
    ]

    for name, _ in scenarios:
        print(f"Generating: {name} ...")

        if name == "stress_fibers":
            actin, gt = make_stress_fibers(cell_mask)
            bact = empty_bact
        elif name == "cortical_ring":
            actin, gt = make_cortical_ring(cell_mask)
            bact = empty_bact
        elif name == "dense_mesh":
            actin, gt = make_dense_mesh(cell_mask)
            bact = empty_bact
        elif name == "infected_cell":
            actin, bact, gt = make_infected_cell(cell_mask, nuc_mask)
        elif name == "empty_control":
            actin, gt = make_empty_control(cell_mask)
            bact = empty_bact

        stack = build_stack(actin, nuc_ch, bf_ch, bact)
        out_path = out_dir / f"{name}.tif"
        tifffile.imwrite(str(out_path), stack, imagej=True,
                         metadata={"axes": "ZCYX"})
        print(f"  Saved: {out_path}  shape={stack.shape}  dtype={stack.dtype}")

        gt["filename"] = f"{name}.tif"
        gt["image_shape"] = str(stack.shape)
        gt["channels"] = "actin=0, nuclei=1, brightfield=2, bacteria=3"
        ground_truth_rows.append(gt)

    # Write ground truth CSV
    gt_path = out_dir / "ground_truth.csv"
    all_keys = []
    for row in ground_truth_rows:
        for k in row:
            if k not in all_keys:
                all_keys.append(k)

    with open(gt_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=all_keys)
        writer.writeheader()
        for row in ground_truth_rows:
            writer.writerow({k: row.get(k, "") for k in all_keys})

    print(f"\nGround truth CSV: {gt_path}")
    print("\nDone. Directory contents:")
    for p in sorted(out_dir.iterdir()):
        size_kb = p.stat().st_size // 1024
        print(f"  {p.name:<30}  {size_kb:>6} KB")

    print("""
Next steps:
  1. Point config_v3_template.json input_dir to the synthetic_test/ folder
  2. Run: python main_pipeline_v3.py
  3. Compare pipeline CSV to ground_truth.csv to measure accuracy

Key accuracy checks:
  stress_fibers  -> mean_tortuosity_edge close to 1.0
  cortical_ring  -> cortical_weighted_density > inner_weighted_density
  infected_cell  -> bact_density_perinuclear > 0, bact_density_cortical ~ 0
  empty_control  -> pipeline completes without crash
""")


if __name__ == "__main__":
    out_dir = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path("synthetic_test")
    main(out_dir)
