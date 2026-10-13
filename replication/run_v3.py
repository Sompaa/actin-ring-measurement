"""Run the old TopoActin V3 batch pipeline on a folder of fields.

    python replication/run_v3.py INPUT_DIR OUTPUT_DIR
    python replication/run_v3.py INPUT_DIR OUTPUT_DIR --measure-on raw
    python replication/run_v3.py INPUT_DIR OUTPUT_DIR --rings   # + 10-ring profile

INPUT_DIR is a folder of *_Z_stack_merged.tif(f) files (see prepare_input.py).
The per-cell CSV is written to OUTPUT_DIR/results/cytoskeleton_analysis_v3_*.csv.

--measure-on adjusted (default) measures intensities on the brightness/contrast-
adjusted actin image, which is what the January 2026 runs behind the paper did.
The code in old_pipeline/ was changed on 10 Feb 2026 to measure on the raw image
instead; that gives intensities 12-25x lower and does not match the paper.
Use --measure-on raw only to see that difference. The old code itself is not
edited: the switches are applied at run time by v3_paper.py.

--rings adds the 10-ring Sato-weighted density profile (ring 0 = nucleus,
rings 1-9 = cytoplasm) and the rings merged 3-3-3 into three zones.
"""
import argparse
import json
from pathlib import Path

from v3_paper import DEFAULT_CONFIG, add_ring_profile, apply_paper_settings, load_config


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_dir")
    ap.add_argument("output_dir")
    ap.add_argument("--config", default=str(DEFAULT_CONFIG))
    ap.add_argument("--measure-on", choices=["adjusted", "raw"], default="adjusted")
    ap.add_argument("--cpu", action="store_true", help="run Cellpose on the CPU (slow)")
    ap.add_argument("--rings", action="store_true",
                    help="also write the 10-ring Sato-weighted density profile")
    args = ap.parse_args()

    cfg = load_config(args.config, cpu=args.cpu)
    cfg["input_dir"] = str(Path(args.input_dir).resolve())
    cfg["output_dir"] = str(Path(args.output_dir).resolve())
    apply_paper_settings(cfg, measure_on=args.measure_on)
    if args.rings:
        add_ring_profile()

    import main_pipeline_v3 as mp

    out = Path(cfg["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    (out / "settings_used.json").write_text(
        json.dumps({"measure_on": args.measure_on, "rings": args.rings, "config": cfg}, indent=1),
        encoding="utf-8")
    print(f"Measuring intensities on the {args.measure_on} image; settings saved to {out / 'settings_used.json'}")
    mp.process_experiment(cfg)


if __name__ == "__main__":
    main()
