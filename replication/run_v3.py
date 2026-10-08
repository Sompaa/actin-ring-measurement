"""Run the old TopoActin V3 batch pipeline on a folder of fields.

    python replication/run_v3.py INPUT_DIR OUTPUT_DIR
    python replication/run_v3.py INPUT_DIR OUTPUT_DIR --measure-on raw

INPUT_DIR is a folder of *_Z_stack_merged.tif(f) files (see prepare_input.py).
The per-cell CSV is written to OUTPUT_DIR/results/cytoskeleton_analysis_v3_*.csv.

--measure-on adjusted (default) measures intensities on the brightness/contrast-
adjusted actin image, which is what the January 2026 runs behind the paper did.
The code in old_pipeline/ was changed on 10 Feb 2026 to measure on the raw image
instead; that gives intensities 12-25x lower and does not match the paper.
Use --measure-on raw only to see that difference. The old code itself is not
edited: the switch is applied here, at run time.
"""
import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "old_pipeline"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_dir")
    ap.add_argument("output_dir")
    ap.add_argument("--config", default=str(REPO / "replication" / "config_paper.json"))
    ap.add_argument("--measure-on", choices=["adjusted", "raw"], default="adjusted")
    ap.add_argument("--cpu", action="store_true", help="run Cellpose on the CPU (slow)")
    args = ap.parse_args()

    import main_pipeline_v3 as mp
    import segmentation_v3

    cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    cfg["input_dir"] = str(Path(args.input_dir).resolve())
    cfg["output_dir"] = str(Path(args.output_dir).resolve())
    if args.cpu:
        cfg.setdefault("processing", {})["use_gpu"] = False

    if args.measure_on == "adjusted":
        original = segmentation_v3.CellSegmenter.segment_image

        def segment_and_measure_adjusted(self, path):
            result = original(self, path)
            if result is not None:
                result["actin_image_raw"] = result["actin_image"]
                if "actin_volume_3d" in result:
                    result["actin_volume_3d_raw"] = result["actin_volume_3d"]
            return result

        segmentation_v3.CellSegmenter.segment_image = segment_and_measure_adjusted
        mp.CellSegmenter = segmentation_v3.CellSegmenter

    # skeletonization_v3 calls skimage's sato() with its default black_ridges=True,
    # which traces the dark gaps between fibres. Actin fibres are bright; the
    # earlier Skeletactin code used black_ridges=False. The config key decides.
    black = cfg.get("skeletonization", {}).get("sato_black_ridges", True)
    if not black:
        import skeletonization_v3
        from skimage import filters as _skf

        def sato_bright(self, volume):
            return _skf.sato(volume, sigmas=self.sato_sigmas, mode="reflect", black_ridges=False)

        skeletonization_v3.CellSkeletonizer._apply_sato_3d = sato_bright
        mp.CellSkeletonizer = skeletonization_v3.CellSkeletonizer

    out = Path(cfg["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    (out / "settings_used.json").write_text(
        json.dumps({"measure_on": args.measure_on, "config": cfg}, indent=1), encoding="utf-8")
    print(f"Measuring intensities on the {args.measure_on} image; settings saved to {out / 'settings_used.json'}")
    mp.process_experiment(cfg)


if __name__ == "__main__":
    main()
