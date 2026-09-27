#!/usr/bin/env python3
"""Export the full CellChatDB.mouse supplementary tables as CSV files."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path("/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis/manuscript/supplementary_package")
MANIFEST = ROOT / "full_cellchat_tables_data.json"
OUT = ROOT / "full_cellchat_tables"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    for table in manifest["tables"]:
        frame = pd.DataFrame(table["rows"], columns=table["columns"])
        filename = table["key"].replace("_", "-") + ".csv"
        frame.to_csv(OUT / filename, index=False, encoding="utf-8-sig")
        print(f"Wrote {filename}: {len(frame)} rows")


if __name__ == "__main__":
    main()
