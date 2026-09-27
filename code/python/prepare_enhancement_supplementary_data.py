#!/usr/bin/env python3
"""Prepare supplementary tables for the optional enhancement analyses."""

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


def combine_enrichment(pattern: str, analysis_name: str) -> pd.DataFrame:
    frames = []
    for gene_set in ("M12_module", "Severity_positive_FDR005"):
        frame = read_csv(ROOT / "enhancement_analysis" / f"{gene_set}_{pattern}.csv")
        frame.insert(0, "gene_set", gene_set)
        frame.insert(1, "analysis_type", analysis_name)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    tables = []

    communication = read_csv(ROOT / "cellchat_enhancement" / "cellchatdb_communication_scores.csv")
    add_table(
        tables,
        "S19_CellChatDB_communication_scores",
        "Supplementary Table S19. CellChatDB-based communication scores.",
        "CellChatDB mouse ligand-receptor communication scores by age, sender cell type, receiver cell type, and interaction.",
        communication,
        "Scores use CellChatDB interaction and complex definitions with a CellChat-style ligand-receptor product score.",
    )

    changes = read_csv(ROOT / "cellchat_enhancement" / "cellchatdb_interaction_changes.csv")
    add_table(
        tables,
        "S20_CellChatDB_interaction_changes",
        "Supplementary Table S20. CellChatDB-based communication changes.",
        "Changes in ligand-receptor communication scores between 24-month and 3-month cochleae.",
        changes,
    )

    tf = combine_enrichment("transcription_factors", "transcription_factors")
    add_table(
        tables,
        "S21_Transcription_factor_enrichment",
        "Supplementary Table S21. Transcription factor enrichment.",
        "ChEA 2022 transcription-factor enrichment for the M12 module and severity-positive gene sets.",
        tf,
    )

    mirna = combine_enrichment("micrornas", "micrornas")
    add_table(
        tables,
        "S22_miRNA_target_enrichment",
        "Supplementary Table S22. miRNA target enrichment.",
        "TargetScan miRNA target enrichment for the M12 module and severity-positive gene sets.",
        mirna,
        "No miRNA term reached FDR < 0.05 in the current analysis.",
    )

    go = combine_enrichment("go_biological_process", "GO Biological Process")
    kegg = combine_enrichment("kegg_mouse", "KEGG mouse")
    functional = pd.concat([go, kegg], ignore_index=True)
    add_table(
        tables,
        "S23_GO_KEGG_enrichment",
        "Supplementary Table S23. Second-pass GO and KEGG enrichment.",
        "GO Biological Process 2025 and KEGG mouse enrichment results for the M12 module and severity-positive gene sets.",
        functional,
    )

    drug = combine_enrichment("drug_repositioning", "drug repositioning")
    add_table(
        tables,
        "S24_Drug_repositioning",
        "Supplementary Table S24. Drug repositioning signatures.",
        "DSigDB enrichment results for the M12 module and severity-positive gene sets.",
        drug,
        "Drug terms are computational perturbation signatures and are not therapeutic recommendations.",
    )

    formal_cellchat = read_csv(ROOT / "cellchat_formal" / "CellChat_all_communications.csv")
    add_table(
        tables,
        "S25_Formal_CellChat_communications",
        "Supplementary Table S25. Formal CellChat communication results.",
        "CellChat 1.6.1 communication probabilities and permutation P values for the severity-associated CellChatDB interaction subset.",
        formal_cellchat,
        "Formal CellChat used CellChatDB.mouse, triMean expression, 100 bootstrap permutations, and a minimum of 10 cells per group.",
    )

    payload = {
        "title": "Supplementary Tables S19-S24",
        "manuscript_title": "Peripheral-Central Convergence of Macrophage and Microglial Programs in Presbycusis Severity",
        "tables": tables,
    }
    output = OUT / "enhancement_tables_data.json"
    output.write_text(json.dumps(payload, allow_nan=False), encoding="utf-8")
    print(f"Wrote {output}")
    for table in tables:
        print(f"{table['key']}: {table['row_count']} rows x {table['column_count']} columns")


if __name__ == "__main__":
    main()
