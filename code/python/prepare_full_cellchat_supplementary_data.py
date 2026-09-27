#!/usr/bin/env python3
"""Prepare the full CellChatDB.mouse supplementary tables (S35-S37)."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path("/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis")
OUT = ROOT / "manuscript" / "supplementary_package"


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, encoding="utf-8-sig")


def clean_records(frame: pd.DataFrame) -> list[dict[str, object]]:
    frame = frame.astype(object).where(pd.notna(frame), None)
    return frame.to_dict(orient="records")


def add_table(tables, key, title, description, frame, notes=""):
    tables.append(
        {
            "key": key,
            "title": title,
            "description": description,
            "notes": notes,
            "rows": clean_records(frame),
            "row_count": int(len(frame)),
            "column_count": int(len(frame.columns)),
            "columns": list(frame.columns),
        }
    )


def build_summary() -> pd.DataFrame:
    frames = []
    for age in (3, 12, 24):
        frame = read_csv(
            ROOT / "cellchat_full" / "final" / f"CellChat_full_summary_{age}M.csv"
        )
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)
    return combined[
        [
            "age_months",
            "n_interactions_evaluated",
            "n_nonzero_edges",
            "n_significant_edges",
            "n_significant_interactions",
            "n_significant_pathways",
            "nboot",
            "seed_use",
        ]
    ]


def main() -> None:
    tables = []

    add_table(
        tables,
        "S35_Full_CellChat_summary",
        "Supplementary Table S35. Full CellChatDB.mouse run summary.",
        "Per-age summary of the complete 2,019-interaction CellChatDB.mouse run in GSE274279, including evaluated interactions, nonzero and significant edges, significant ligand-receptor interactions, significant pathways, and permutation settings.",
        build_summary(),
        "All 2,019 mouse ligand-receptor pairs were evaluated per age. Genes absent from the GSE274279 matrix were represented as all-zero rows so that every database interaction could be tested; such interactions return zero probability.",
    )
    add_table(
        tables,
        "S36_Full_CellChat_significant_interactions",
        "Supplementary Table S36. Significant ligand-receptor interactions in the full CellChatDB.mouse run.",
        "Cell-type-resolved significant (permutation p < 0.05) ligand-receptor edges across the three ages.",
        read_csv(
            ROOT / "cellchat_full" / "final" / "CellChat_full_significant_all_ages.csv"
        ),
        "P values are CellChat permutation values computed from 100 bootstrap replicates and are reported without additional multiple-testing correction.",
    )
    add_table(
        tables,
        "S37_Full_CellChat_pathways",
        "Supplementary Table S37. Pathway-level communication in the full CellChatDB.mouse run.",
        "Cell-type-resolved pathway-level communication probabilities aggregated from significant ligand-receptor edges across the three ages.",
        read_csv(
            ROOT / "cellchat_full" / "final" / "CellChat_full_pathway_all_ages.csv"
        ),
    )

    payload = {
        "title": "Supplementary Tables S35-S37",
        "manuscript_title": "Peripheral-Central Convergence of Macrophage and Microglial Programs in Presbycusis Severity",
        "tables": tables,
    }
    output = OUT / "full_cellchat_tables_data.json"
    output.write_text(json.dumps(payload, allow_nan=False), encoding="utf-8")
    print(f"Wrote {output}")
    for table in tables:
        print(f"{table['key']}: {table['row_count']} rows x {table['column_count']} columns")


if __name__ == "__main__":
    main()
