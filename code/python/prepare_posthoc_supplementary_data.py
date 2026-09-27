#!/usr/bin/env python3
"""Prepare sensitivity and myeloid-subclustering supplementary tables."""

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


def build_model_sensitivity() -> pd.DataFrame:
    frames = []
    for dataset, model, filename in [
        ("Cochlea", "Severity + sex", "cochlea_severity_sex_adjusted.csv"),
        ("Cochlea", "Severity + age group", "cochlea_severity_age_adjusted.csv"),
        ("Inferior colliculus", "Severity + sex", "central_severity_sex_adjusted.csv"),
        ("Inferior colliculus", "Severity + age group", "central_severity_age_adjusted.csv"),
    ]:
        frame = read_csv(ROOT / "sensitivity_analysis" / filename)
        frame.insert(0, "model", model)
        frame.insert(0, "dataset", dataset)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    tables = []

    add_table(
        tables,
        "S26_Model_sensitivity",
        "Supplementary Table S26. Severity-model sensitivity results.",
        "Complete limma results after adjustment for sex or age group in GSE49543 and GSE49522.",
        build_model_sensitivity(),
    )
    add_table(
        tables,
        "S27_Candidate_model_sensitivity",
        "Supplementary Table S27. Candidate-gene model sensitivity.",
        "Candidate gene coefficients from severity-only, sex-adjusted, and age-adjusted models.",
        read_csv(ROOT / "sensitivity_analysis" / "candidate_models_sensitivity.csv"),
    )
    add_table(
        tables,
        "S28_Candidate_bootstrap",
        "Supplementary Table S28. Candidate-gene bootstrap estimates.",
        "Bootstrap confidence intervals for candidate-gene severity coefficients.",
        read_csv(ROOT / "sensitivity_analysis" / "candidate_severity_bootstrap.csv"),
    )
    add_table(
        tables,
        "S29_Myeloid_cell_metadata",
        "Supplementary Table S29. Myeloid cell metadata.",
        "Cell-level age, cluster, and myeloid state assignments for GSE274279.",
        read_csv(ROOT / "myeloid_subclustering" / "myeloid_cell_metadata.csv"),
    )
    add_table(
        tables,
        "S30_Myeloid_cluster_markers",
        "Supplementary Table S30. Myeloid cluster markers.",
        "Positive markers for each myeloid subcluster.",
        read_csv(ROOT / "myeloid_subclustering" / "myeloid_cluster_markers.csv"),
    )
    add_table(
        tables,
        "S31_Myeloid_module_scores",
        "Supplementary Table S31. Myeloid module scores by state and age.",
        "Mean homeostatic, MHC-II, complement, phagocytic, interferon, APP-CD74, and PTPRC-MRC1 scores.",
        read_csv(ROOT / "myeloid_subclustering" / "myeloid_module_scores_by_state_age.csv"),
    )
    add_table(
        tables,
        "S32_Myeloid_bootstrap",
        "Supplementary Table S32. Myeloid module-score bootstrap.",
        "Equal-cell bootstrap estimates of age-related module-score changes.",
        read_csv(ROOT / "myeloid_subclustering" / "myeloid_module_bootstrap_samples.csv"),
    )
    add_table(
        tables,
        "S33_Myeloid_cluster_proportions",
        "Supplementary Table S33. Myeloid cluster proportions.",
        "Cluster proportions by age.",
        read_csv(ROOT / "myeloid_subclustering" / "myeloid_cluster_proportions_by_age.csv"),
    )
    add_table(
        tables,
        "S34_Myeloid_UMAP",
        "Supplementary Table S34. Myeloid UMAP coordinates.",
        "UMAP coordinates, state assignment, and age for 147 myeloid nuclei.",
        read_csv(ROOT / "myeloid_subclustering" / "myeloid_umap_coordinates.csv"),
    )

    payload = {
        "title": "Supplementary Tables S26-S34",
        "manuscript_title": "Peripheral-Central Convergence of Macrophage and Microglial Programs in Presbycusis Severity",
        "tables": tables,
    }
    output = OUT / "posthoc_tables_data.json"
    output.write_text(json.dumps(payload, allow_nan=False), encoding="utf-8")
    print(f"Wrote {output}")
    for table in tables:
        print(f"{table['key']}: {table['row_count']} rows x {table['column_count']} columns")


if __name__ == "__main__":
    main()
