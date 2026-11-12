#!/usr/bin/env python3
"""Add display columns for exact-zero P values without changing numeric fields."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


P_THRESHOLD = "<2.2e-16"
P_COLUMNS = ("p_value", "pval")
FDR_COLUMNS = ("fdr_bh",)


def display_value(value: str) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return value
    return P_THRESHOLD if numeric == 0 else value


def add_display_columns(path: Path) -> tuple[int, int]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        rows = list(reader)

    additions: list[tuple[str, str]] = []
    for column in P_COLUMNS:
        if column in fieldnames and f"{column}_display" not in fieldnames:
            additions.append((column, f"{column}_display"))
    for column in FDR_COLUMNS:
        if column in fieldnames and f"{column}_display" not in fieldnames:
            additions.append((column, f"{column}_display"))

    if not additions:
        return 0, len(rows)

    output_fields = fieldnames + [display for _, display in additions]
    for row in rows:
        for source, display in additions:
            row[display] = display_value(row.get(source, ""))

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields)
        writer.writeheader()
        writer.writerows(rows)
    return len(additions), len(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", type=Path)
    args = parser.parse_args()

    for path in args.files:
        columns, rows = add_display_columns(path)
        print(f"{path}: added_display_columns={columns} rows={rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
