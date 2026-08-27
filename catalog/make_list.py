#!/usr/bin/env python3

import json
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--dataset", type=str, required=True)
parser.add_argument("-s", "--sample", type=str, required=True)
parser.add_argument("-v", "--version", type=int, default=-1)
parser.add_argument("-t", "--date", type=str, default=None,
                     help="Select entry by date (e.g. 08_21_2026) instead of version")
parser.add_argument("-n", "--number", type=int, default=100)

args = parser.parse_args()

with open("catalog.json") as f:
    catalog = json.load(f)

try:
    entries = catalog[args.dataset][args.sample]
except KeyError:
    raise SystemExit(
        f"Unknown dataset/sample '{args.dataset}/{args.sample}'"
    )

if args.date is not None:
    dates = sorted(x["date"] for x in entries)
    try:
        entry = next(x for x in entries if x["date"] == args.date)
    except StopIteration:
        raise SystemExit(
            f"Date {args.date} not found.\n"
            f"Available dates: {', '.join(dates)}"
        )
else:
    versions = sorted(x["version"] for x in entries)

    # Default to latest
    if args.version == -1:
        args.version = versions[-1]

    try:
        entry = next(x for x in entries if x["version"] == args.version)
    except StopIteration:
        raise SystemExit(
            f"Version {args.version} not found.\n"
            f"Available versions: {', '.join(map(str, versions))}"
        )

files = sorted(Path(entry["path"]).glob("*.root"))

outdir = Path(args.sample)
outdir.mkdir(exist_ok=True)

for stale in outdir.glob(f"{args.sample}-*.list"):
    stale.unlink()

for i in range(0, len(files), args.number):
    outfile = outdir / f"{args.sample}-{i//args.number:03d}.list"

    with open(outfile, "w") as f:
        for file in files[i:i+args.number]:
            f.write(f"{file}\n")

print(
    f"Wrote {len(files)} files from version {entry['version']} ({entry['date']})"
)