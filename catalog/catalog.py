#!/usr/bin/env python3

import json
import re
from pathlib import Path

# load config file
config_file = Path(__file__).parent / "config"
with open(config_file) as f:
    config = {}
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        key, value = line.split("=", 1)
        config[key.strip()] = value.strip()
if "STORAGE_DIR" not in config or "CATALOG_FILE" not in config:
    raise ValueError("STORAGE_DIR and CATALOG_FILE must be defined in config file")

ROOT = Path(config["STORAGE_DIR"])
OUT  = Path(config["CATALOG_FILE"])
print(f"Using storage directory: {ROOT}")
print(f"Using catalog file: {OUT}")

# production dirs are named e.g. "08_14_2026_v001" -> date=08_14_2026, version=1
VERSION_DIR_RE = re.compile(r"^(?P<date>\d{2}_\d{2}_\d{4})_v(?P<version>\d+)$")


def count_root_files(path):
    return sum(1 for _ in path.glob("*.root"))


catalog = {}

for version_dir in sorted(p for p in ROOT.iterdir() if p.is_dir()):
    match = VERSION_DIR_RE.match(version_dir.name)
    if not match:
        print(f"Skipping unrecognized directory: {version_dir.name}")
        continue

    date = match.group("date")
    version = int(match.group("version"))

    # trees/: flat directory of *.root files, plus optional
    # trees/<sample>/{scaled,unsclaed}/*.root (pass3 output, mirrors dsts/)
    trees_dir = version_dir / "trees"
    if trees_dir.is_dir():
        catalog.setdefault("trees", {}).setdefault("trees", []).append({
            "version": version,
            "date": date,
            "path": str(trees_dir),
            "nfiles": count_root_files(trees_dir),
        })

        for sample_dir in sorted(p for p in trees_dir.iterdir() if p.is_dir()):
            for kind_dir in sorted(p for p in sample_dir.iterdir() if p.is_dir()):
                sample_name = f"{sample_dir.name}_{kind_dir.name}"
                catalog.setdefault("trees", {}).setdefault(sample_name, []).append({
                    "version": version,
                    "date": date,
                    "path": str(kind_dir),
                    "nfiles": count_root_files(kind_dir),
                })

    # dsts/<sample>/{scaled,unsclaed}/*.root
    dsts_dir = version_dir / "dsts"
    if dsts_dir.is_dir():
        for sample_dir in sorted(p for p in dsts_dir.iterdir() if p.is_dir()):
            for kind_dir in sorted(p for p in sample_dir.iterdir() if p.is_dir()):
                sample_name = f"{sample_dir.name}_{kind_dir.name}"
                catalog.setdefault("dsts", {}).setdefault(sample_name, []).append({
                    "version": version,
                    "date": date,
                    "path": str(kind_dir),
                    "nfiles": count_root_files(kind_dir),
                })

for dataset in catalog.values():
    for versions in dataset.values():
        versions.sort(key=lambda x: x["version"])

with open(OUT, "w") as f:
    json.dump(catalog, f, indent=2)

print(f"Wrote {OUT}")
