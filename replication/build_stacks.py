"""Rebuild Z_Stack_Merged/ from a raw Harmony (Opera Phenix) export.

    python replication/build_stacks.py PLATE_DIR
    python replication/build_stacks.py PLATE_DIR --fields replication/field_lists/Exp10_representative_fields.txt

PLATE_DIR is the measurement folder holding Images/Index.xml and the raw
single-plane TIFFs (r03c09f13p01-ch3sk1fk1fl1.tiff ...). For every field it
writes PLATE_DIR/Z_Stack_Merged/<Well>_F<nn>_Z_stack_merged.tiff: an ImageJ
hyperstack (Z, C, Y, X), uint16, with the actin channel first and the nucleus
channel second. Pixel values are copied unchanged. This reproduces the stacks
the paper runs used (made in ImageJ 1.54p), which were checked to hold the raw
planes byte for byte.

Channel order matters: the V3 pipeline reads channel 0 as actin and channel 1
as nuclei. Harmony numbers the channels differently per plate, so they are
chosen by name from Index.xml (override with --actin-ch / --nuclei-ch).
Brightfield is left out, as in the original stacks.

Existing stacks are not overwritten unless --overwrite is given.
"""
import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import tifffile

ACTIN_NAMES = ("cellmask", "phalloidin", "actin", "alexa 488", "alexa_488")
NUCLEI_NAMES = ("hoechst", "dapi")
RAW = re.compile(r"r(\d{2})c(\d{2})f(\d{2,3})p(\d{2})-ch(\d+)", re.I)


def channels_from_index(index_xml):
    text = Path(index_xml).read_text(encoding="utf-8", errors="ignore")
    found = {}
    for block in re.findall(r"<Entry ChannelID=\"(\d+)\"[^>]*>(.*?)</Entry>", text, re.S):
        cid, body = block
        m = re.search(r"<ChannelName>([^<]*)</ChannelName>", body)
        if m:
            found[int(cid)] = m.group(1)
    return found


def pick(found, names):
    for cid, name in sorted(found.items()):
        if any(n in name.lower() for n in names):
            return cid
    return None


def read_fields(path):
    out = set()
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.add(line)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plate_dir")
    ap.add_argument("--fields", help="only build the fields in this list file")
    ap.add_argument("--actin-ch", type=int)
    ap.add_argument("--nuclei-ch", type=int)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    plate = Path(args.plate_dir)
    images = plate / "Images"
    index = next((p for p in (images / "Index.xml", images / "index.xml") if p.exists()), None)
    found = channels_from_index(index) if index else {}
    actin = args.actin_ch or pick(found, ACTIN_NAMES)
    nuclei = args.nuclei_ch or pick(found, NUCLEI_NAMES)
    if not actin or not nuclei:
        sys.exit(f"Could not tell the channels apart from Index.xml ({found}); pass --actin-ch and --nuclei-ch")
    print(f"actin = ch{actin} ({found.get(actin, '?')}), nuclei = ch{nuclei} ({found.get(nuclei, '?')})")

    planes = defaultdict(dict)
    for f in images.glob("*.tif*"):
        m = RAW.match(f.name)
        if m:
            row, col, fld, z, ch = (int(g) for g in m.groups())
            name = f"{chr(ord('A') + row - 1)}{col}_F{fld:02d}"
            planes[name][(z, ch)] = f
    wanted = read_fields(args.fields) if args.fields else set(planes)

    out_dir = plate / "Z_Stack_Merged"
    out_dir.mkdir(exist_ok=True)
    made = skipped = 0
    for name in sorted(wanted):
        if name not in planes:
            print(f"  {name}: no raw planes found, skipped")
            continue
        dst = out_dir / f"{name}_Z_stack_merged.tiff"
        if dst.exists() and not args.overwrite:
            skipped += 1
            continue
        zs = sorted({z for z, _ in planes[name]})
        missing = [(z, c) for z in zs for c in (actin, nuclei) if (z, c) not in planes[name]]
        if missing:
            print(f"  {name}: missing planes {missing}, skipped")
            continue
        stack = np.stack([np.stack([tifffile.imread(planes[name][(z, c)]) for c in (actin, nuclei)])
                          for z in zs]).astype(np.uint16)
        tifffile.imwrite(dst, stack, imagej=True, metadata={"axes": "ZCYX", "mode": "composite"})
        made += 1
    print(f"{made} stack(s) written to {out_dir}, {skipped} already there (use --overwrite to rebuild)")


if __name__ == "__main__":
    main()
