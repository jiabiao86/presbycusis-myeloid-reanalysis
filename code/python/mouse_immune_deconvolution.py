#!/usr/bin/env python3
from __future__ import annotations

import gzip
import json
import math
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = Path("/Users/jijiabiao/Desktop/耳聋课题")
OUT = ROOT / "2026_supplement_analysis"
TEMP = Path("/tmp/ear_study_learning")
EXPRESSION = ROOT / "2022耳聋耳鸣" / "public" / "exp.txt"

MARKERS = {
    "Macrophage/Microglia": [
        "Cd68",
        "C1qa",
        "C1qb",
        "C1qc",
        "Tyrobp",
        "Csf1r",
        "Fcgr3",
        "Mpeg1",
        "Ms4a7",
        "Laptm5",
        "Ctss",
        "Apoe",
        "Cd74",
        "H2-Aa",
        "H2-Eb1",
    ],
    "Microglia-like": ["P2ry12", "Tmem119", "Cx3cr1", "Hexb", "Trem2"],
    "Monocytes": ["Ly6c2", "Ccr2", "F13a1", "Plac8", "Ms4a6c"],
    "Neutrophils": ["S100a8", "S100a9", "Retnlg", "Il1b", "Mpo"],
    "T cells": ["Cd3d", "Cd3e", "Cd3g", "Lck", "Ms4a4b"],
    "B cells": ["Cd79a", "Cd79b", "Ms4a1"],
    "NK cells": ["Nkg7", "Klrb1c", "Ncr1"],
    "Dendritic cells": ["Xcr1", "Itgax", "Cd209a", "H2-Ab1"],
    "Mast cells": ["Cpa3", "Kit", "Tpsab1"],
}


def bh(p_values: np.ndarray) -> np.ndarray:
    p_values = np.clip(np.asarray(p_values, dtype=float), np.nextafter(0.0, 1.0), 1.0)
    count = len(p_values)
    order = np.argsort(p_values)
    adjusted = p_values[order] * count / (np.arange(count) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(count, dtype=float)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def marker_scores(
    expression: pd.DataFrame,
    markers: dict[str, list[str]],
) -> pd.DataFrame:
    rows = {}
    for cell_type, genes in markers.items():
        available = [gene for gene in genes if gene in expression.index]
        if len(available) < 2:
            continue
        sub = expression.loc[available]
        standardized = sub.sub(sub.mean(axis=1), axis=0).div(
            sub.std(axis=1).replace(0, np.nan),
            axis=0,
        )
        rows[cell_type] = standardized.mean(axis=0)
    return pd.DataFrame(rows, index=expression.columns)


def correlate_with_severity(
    scores: pd.DataFrame,
    groups: list[str],
    dataset: str,
    gene_universe: pd.Index,
) -> list[dict]:
    severity = np.array([{"YC": 0, "MA": 1, "MP": 2, "SP": 3}[group] for group in groups], dtype=float)
    rows = []
    for cell_type in scores.columns:
        x = severity - severity.mean()
        y = scores[cell_type].to_numpy(dtype=float)
        y = y - y.mean()
        denominator = math.sqrt(float((x**2).sum()) * float((y**2).sum()))
        r = 0.0 if denominator == 0 else float((x * y).sum() / denominator)
        if abs(r) >= 1:
            p_value = 0.0
        else:
            t_value = r * math.sqrt((len(groups) - 2) / max(1 - r**2, 1e-300))
            p_value = math.erfc(abs(t_value) / math.sqrt(2))
        rows.append(
            {
                "dataset": dataset,
                "cell_type": cell_type,
                "correlation_with_severity": r,
                "p_value": p_value,
                "genes_used": len(
                    [gene for gene in MARKERS[cell_type] if gene in gene_universe]
                ),
            }
        )
    result = pd.DataFrame(rows)
    result["fdr_bh"] = bh(result["p_value"].to_numpy())
    return result.to_dict(orient="records")


def load_gse49543() -> tuple[pd.DataFrame, list[str]]:
    expression = pd.read_csv(EXPRESSION, sep="\t", index_col=0)
    expression = expression.groupby(level=0).mean()
    expression = np.log2(expression + 1.0)
    groups = [column.split("-")[0] for column in expression.columns]
    return expression, groups


def load_gse49522() -> tuple[pd.DataFrame, list[str]]:
    lines = []
    status_line = None
    in_table = False
    with gzip.open(TEMP / "GSE49522_series_matrix.txt.gz", "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!Sample_characteristics_ch1") and "hearing status:" in line:
                status_line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                lines.append(line.rstrip("\n"))
    if status_line is None:
        raise RuntimeError("Missing GSE49522 status")
    statuses = [value.strip().strip('"').split(":", 1)[-1].strip() for value in status_line.split("\t")[1:]]
    mapping = {
        "Young Control": "YC",
        "Middle-Aged": "MA",
        "Mild Presbycusis": "MP",
        "Severe Presbycusis": "SP",
    }
    groups = [mapping[status] for status in statuses]
    platform_rows = []
    in_platform = False
    with (TEMP / "GPL339_full.txt").open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_platform = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_platform:
                platform_rows.append(line.rstrip("\n"))
    platform = pd.read_csv(StringIO("\n".join(platform_rows)), sep="\t", dtype=str)
    mapping_frame = platform[["ID", "Gene Symbol"]].dropna()
    mapping_frame = mapping_frame[
        (mapping_frame["Gene Symbol"] != "")
        & (mapping_frame["Gene Symbol"] != "---")
        & (~mapping_frame["Gene Symbol"].str.contains("///", regex=False))
    ]
    probe_to_symbol = dict(zip(mapping_frame["ID"], mapping_frame["Gene Symbol"]))
    matrix = pd.read_csv(StringIO("\n".join(lines)), sep="\t", index_col=0)
    matrix.index = [probe_to_symbol.get(str(index), "") for index in matrix.index]
    matrix = matrix[matrix.index != ""].groupby(level=0).mean()
    return np.log2(matrix + 1.0), groups


def load_gse233798() -> tuple[pd.DataFrame, np.ndarray]:
    frame = pd.read_excel(TEMP / "GSE233798_FPKM.xlsx")
    frame = frame.dropna(subset=["gene_symbol"]).copy()
    frame["gene_symbol"] = frame["gene_symbol"].astype(str)
    columns = [column for column in frame.columns if column.startswith("fpkm_")]
    frame[columns] = frame[columns].apply(pd.to_numeric, errors="coerce")
    frame["mean_expression"] = frame[columns].mean(axis=1)
    frame = (
        frame.sort_values("mean_expression", ascending=False)
        .drop_duplicates("gene_symbol")
        .set_index("gene_symbol")
    )
    young = [column for column in columns if "young" in column.lower()]
    old = [column for column in columns if "old" in column.lower()]
    return np.log2(frame[old + young] + 1.0), np.array([1] * len(old) + [0] * len(young))


def main() -> int:
    cochlea, cochlea_groups = load_gse49543()
    central, central_groups = load_gse49522()
    external, external_old = load_gse233798()

    cochlea_scores = marker_scores(cochlea, MARKERS)
    central_scores = marker_scores(central, MARKERS)
    external_scores = marker_scores(external, MARKERS)
    cochlea_scores.to_csv(OUT / "mouse_immune_scores_GSE49543.csv", encoding="utf-8-sig")
    central_scores.to_csv(OUT / "mouse_immune_scores_GSE49522.csv", encoding="utf-8-sig")
    external_scores.to_csv(OUT / "mouse_immune_scores_GSE233798.csv", encoding="utf-8-sig")

    rows = []
    rows.extend(
        correlate_with_severity(
            cochlea_scores,
            cochlea_groups,
            "GSE49543",
            cochlea.index,
        )
    )
    rows.extend(
        correlate_with_severity(
            central_scores,
            central_groups,
            "GSE49522",
            central.index,
        )
    )

    external_rows = []
    for cell_type in external_scores.columns:
        young = external_scores.loc[external_old == 0, cell_type].to_numpy(dtype=float)
        old = external_scores.loc[external_old == 1, cell_type].to_numpy(dtype=float)
        statistic, p_value = mannwhitneyu(old, young, alternative="two-sided")
        external_rows.append(
            {
                "dataset": "GSE233798",
                "cell_type": cell_type,
                "old_mean": float(old.mean()),
                "young_mean": float(young.mean()),
                "difference": float(old.mean() - young.mean()),
                "p_value": float(p_value),
                "mannwhitney_u": float(statistic),
            }
        )
    external_result = pd.DataFrame(external_rows)
    external_result["fdr_bh"] = bh(external_result["p_value"].to_numpy())

    summary = pd.DataFrame(rows)
    summary["fdr_bh"] = bh(summary["p_value"].to_numpy())
    summary.to_csv(
        OUT / "mouse_immune_deconvolution_severity.csv",
        index=False,
        encoding="utf-8-sig",
    )
    external_result.to_csv(
        OUT / "mouse_immune_deconvolution_external.csv",
        index=False,
        encoding="utf-8-sig",
    )

    old_scores = pd.read_csv(
        ROOT
        / "figures-中间图形"
        / "figures-中间图形"
        / "fig2"
        / "F"
        / "CIBERSORT_cell_type_FPKM_group.txt",
        sep="\t",
        index_col=0,
    )
    old_scores = old_scores.T
    overlap = old_scores.join(
        cochlea_scores,
        how="inner",
        lsuffix="_CIBERSORT",
        rsuffix="_marker",
    )
    correlation = overlap.corr(method="pearson")
    correlation.to_csv(
        OUT / "mouse_marker_CIBERSORT_correlation.csv",
        encoding="utf-8-sig",
    )

    print(summary.to_string(index=False))
    print("\n", external_result.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
