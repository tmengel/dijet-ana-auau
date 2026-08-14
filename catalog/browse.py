#!/usr/bin/env python3

import json
from pathlib import Path

config_file= Path(__file__).parent / "config"
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

with open(OUT) as f:
    catalog = json.load(f)

for dataset in catalog:

    print(dataset)

    for sample in sorted(catalog[dataset]):

        print(f"  {sample}")

        for v in catalog[dataset][sample]:
            print(
                f"    v{v['version']:02d}"
                f"  {v['date']}"
                f"  ({v['nfiles']} files)"
            )