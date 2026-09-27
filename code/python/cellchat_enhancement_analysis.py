#!/usr/bin/env python3
"""CellChatDB-based ligand-receptor communication scoring for GSE274279."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse


ROOT = Path("/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis")
OUT = ROOT / "cellchat_enhancement"
DATA = Path("/tmp/ear_study_learning/GSE274279") if Path("/tmp/ear_study_learning/GSE274279").exists() else ROOT / "GSE274279"
UMAP_FILE = ROOT / "GSE274279_umap_coordinates.csv"
INTERACTION_FILE = OUT / "cellchatdb_mouse_interactions.csv"
COMPLEX_FILE = OUT / "cellchatdb_mouse_complexes.csv"

AGE_FILES = {
    3: ("3M_matrix.mtx.gz", "3M_barcodes.tsv.gz", "-0"),
    12: ("12M_matrix.mtx.gz", "12M_barcodes.tsv.gz", "-1"),
    24: ("24M_matrix.mtx.gz", "24M_barcodes.tsv.gz", "-2"),
}
CELL_ORDER = [
    "Supporting cells",
    "Fibrocytes",
    "Glia/Schwann",
    "Outer hair cells",
    "Inner hair cells",
    "Macrophages/Microglia",
    "Spiral ganglion neurons",
]
K_HILL = 0.5


def read_features() -> pd.Series:
    features = pd.read_csv(DATA / "features.tsv.gz", sep="\t", header=None, names=["gene_id", "symbol", "type"])
    return features["symbol"].astype(str)


def read_complexes() -> dict[str, list[str]]:
    frame = pd.read_csv(COMPLEX_FILE, encoding="utf-8-sig")
    result: dict[str, list[str]] = {}
    for _, row in frame.iterrows():
        subunits = [
            str(row[column]).strip()
            for column in ("subunit_1", "subunit_2", "subunit_3", "subunit_4")
            if column in row and pd.notna(row[column]) and str(row[column]).strip()
        ]
        result[str(row["complex_name"])] = subunits
    return result


def entity_subunits(entity: str, complex_map: dict[str, list[str]]) -> list[str]:
    entity = str(entity)
    if entity in complex_map:
        return complex_map[entity]
    return [entity]


def aggregate_expression(matrix: sparse.spmatrix, gene_symbols: pd.Series, cell_types: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    present_types = [cell_type for cell_type in CELL_ORDER if cell_type in set(cell_types)]
    means = np.zeros((matrix.shape[0], len(present_types)), dtype=float)
    percents = np.zeros_like(means)
    labels = cell_types.to_numpy()
    for column_index, cell_type in enumerate(present_types):
        mask = labels == cell_type
        subset = matrix[:, mask]
        means[:, column_index] = np.asarray(subset.mean(axis=1)).ravel()
        percents[:, column_index] = np.asarray((subset > 0).sum(axis=1)).ravel() / max(mask.sum(), 1)
    mean_frame = pd.DataFrame(means, index=gene_symbols.to_numpy(), columns=present_types)
    percent_frame = pd.DataFrame(percents * 100.0, index=gene_symbols.to_numpy(), columns=present_types)
    mean_frame = mean_frame.groupby(level=0).mean()
    percent_frame = percent_frame.groupby(level=0).mean()
    return mean_frame, percent_frame


def normalized_expression(path: Path, barcodes: list[str], gene_symbols: pd.Series, required_genes: set[str]) -> tuple[sparse.csr_matrix, pd.Series]:
    matrix = mmread(gzip.open(path, "rb")).tocsr()
    matrix = matrix.astype(float)
    matrix = matrix[gene_symbols.index, :]
    totals = np.asarray(matrix.sum(axis=0)).ravel()
    totals[totals <= 0] = 1.0
    matrix = matrix.multiply((1e4 / totals)[None, :]).tocsr()
    matrix.data = np.log1p(matrix.data)
    keep = gene_symbols.isin(required_genes).to_numpy()
    matrix = matrix[keep, :]
    return matrix, gene_symbols[keep].reset_index(drop=True)


def communication_scores(mean_frame: pd.DataFrame, percent_frame: pd.DataFrame, interactions: pd.DataFrame, complex_map: dict[str, list[str]], age: int) -> pd.DataFrame:
    rows = []
    available = set(mean_frame.index)
    for _, interaction in interactions.iterrows():
        ligands = entity_subunits(interaction["ligand"], complex_map)
        receptors = entity_subunits(interaction["receptor"], complex_map)
        if not all(gene in available for gene in ligands + receptors):
            continue
        ligand_values = mean_frame.loc[ligands].min(axis=0)
        receptor_values = mean_frame.loc[receptors].min(axis=0)
        ligand_percent = percent_frame.loc[ligands].min(axis=0)
        receptor_percent = percent_frame.loc[receptors].min(axis=0)
        for sender in mean_frame.columns:
            for receiver in mean_frame.columns:
                score = (
                    np.sqrt(max(float(ligand_values[sender]), 0.0) * max(float(receptor_values[receiver]), 0.0))
                    * np.sqrt(max(float(ligand_percent[sender]), 0.0) / 100.0 * max(float(receptor_percent[receiver]), 0.0) / 100.0)
                )
                if score <= 0:
                    continue
                rows.append(
                    {
                        "age_months": age,
                        "sender": sender,
                        "receiver": receiver,
                        "interaction_name": interaction["interaction_name"],
                        "pathway_name": interaction["pathway_name"],
                        "ligand": interaction["ligand"],
                        "receptor": interaction["receptor"],
                        "annotation": interaction["annotation"],
                        "ligand_mean": float(ligand_values[sender]),
                        "receptor_mean": float(receptor_values[receiver]),
                        "ligand_percent": float(ligand_percent[sender]),
                        "receptor_percent": float(receptor_percent[receiver]),
                        "communication_score": float(score),
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    interactions = pd.read_csv(INTERACTION_FILE, encoding="utf-8-sig")
    complexes = read_complexes()
    umap = pd.read_csv(UMAP_FILE, index_col=0)
    feature_symbols = read_features()

    required_genes: set[str] = set()
    for _, interaction in interactions.iterrows():
        required_genes.update(entity_subunits(interaction["ligand"], complexes))
        required_genes.update(entity_subunits(interaction["receptor"], complexes))

    age_frames = []
    for age, (matrix_name, barcode_name, suffix) in AGE_FILES.items():
        with gzip.open(DATA / barcode_name, "rt") as handle:
            barcodes = [line.strip() for line in handle if line.strip()]
        index = [f"{barcode}{suffix}" for barcode in barcodes]
        age_umap = umap[umap["age_months"] == age].copy()
        selected = [barcode for barcode in index if barcode in age_umap.index]
        column_index = [index.index(barcode) for barcode in selected]
        matrix, symbols = normalized_expression(DATA / matrix_name, selected, feature_symbols, required_genes)
        matrix = matrix[:, column_index]
        mean_frame, percent_frame = aggregate_expression(matrix, symbols, age_umap.loc[selected, "cell_type"])
        age_scores = communication_scores(mean_frame, percent_frame, interactions, complexes, age)
        age_frames.append(age_scores)
        print(f"Age {age}M: {len(selected)} cells, {len(age_scores)} interaction scores")

    all_scores = pd.concat(age_frames, ignore_index=True)
    all_scores.to_csv(OUT / "cellchatdb_communication_scores.csv", index=False, encoding="utf-8-sig")

    pathway_scores = (
        all_scores.groupby(["age_months", "sender", "receiver", "pathway_name"], as_index=False)["communication_score"]
        .sum()
    )
    pathway_scores.to_csv(OUT / "cellchatdb_pathway_scores.csv", index=False, encoding="utf-8-sig")

    comparison = all_scores.pivot_table(
        index=["sender", "receiver", "interaction_name", "pathway_name", "ligand", "receptor", "annotation"],
        columns="age_months",
        values="communication_score",
        aggfunc="mean",
    ).reset_index()
    comparison.columns.name = None
    comparison["delta_24m_minus_3m"] = comparison.get(24, 0) - comparison.get(3, 0)
    comparison["log2FC_24m_vs_3m"] = np.log2((comparison.get(24, 0) + 1e-6) / (comparison.get(3, 0) + 1e-6))
    comparison["involving_myeloid"] = comparison["sender"].eq("Macrophages/Microglia") | comparison["receiver"].eq("Macrophages/Microglia")
    comparison.sort_values("delta_24m_minus_3m", ascending=False).to_csv(
        OUT / "cellchatdb_interaction_changes.csv",
        index=False,
        encoding="utf-8-sig",
    )

    myeloid = comparison[comparison["involving_myeloid"]].sort_values("delta_24m_minus_3m", ascending=False)
    myeloid.head(100).to_csv(
        OUT / "cellchatdb_top_myeloid_changes.csv",
        index=False,
        encoding="utf-8-sig",
    )
    print(f"Wrote communication tables to {OUT}")


if __name__ == "__main__":
    main()
