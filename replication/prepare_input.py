"""Collect the fields to analyse from a plate into one input folder.

The old pipeline analyses every TIFF it finds under its input folder, so the
way to choose fields is to give it a folder holding only those fields.

    python replication/prepare_input.py PLATE_DIR FIELD_LIST OUT_DIR
    python replication/prepare_input.py PLATE_DIR all OUT_DIR

PLATE_DIR   the plate folder (the one holding Z_Stack_Merged/) or Z_Stack_Merged itself
FIELD_LIST  a file from replication/field_lists/ (one well_field per line, e.g. C9_F13),
            or "all" for every field on the plate
OUT_DIR     created if missing; gets one link (or copy, with --copy) per field

Links are used by default so nothing is duplicated. Use --copy on Windows or
when the plate sits on a drive that will not stay mounted.
"""
import argparse
import os
import shutil
import sys
from pathlib import Path


def read_fields(path):
    names = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            names.append(line)
    return names


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plate_dir")
    ap.add_argument("field_list")
    ap.add_argument("out_dir")
    ap.add_argument("--copy", action="store_true", help="copy the TIFFs instead of linking them")
    args = ap.parse_args()

    plate = Path(args.plate_dir)
    stacks = plate / "Z_Stack_Merged" if (plate / "Z_Stack_Merged").is_dir() else plate
    available = {}
    for f in sorted(stacks.glob("*.tif*")):
        available[f.name.split("_Z_")[0]] = f
    if not available:
        sys.exit(f"No *_Z_stack_merged.tif(f) files found in {stacks}")

    wanted = list(available) if args.field_list == "all" else read_fields(args.field_list)
    missing = [w for w in wanted if w not in available]
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    n = 0
    for name in wanted:
        if name not in available:
            continue
        src, dst = available[name], out / available[name].name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        if args.copy:
            shutil.copy2(src, dst)
        else:
            os.symlink(src.resolve(), dst)
        n += 1

    print(f"{n} field(s) placed in {out}  (plate has {len(available)})")
    if missing:
        print(f"WARNING: {len(missing)} listed field(s) are not on this plate: {', '.join(missing[:10])}"
              + (" ..." if len(missing) > 10 else ""))


if __name__ == "__main__":
    main()
