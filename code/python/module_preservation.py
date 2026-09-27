#!/usr/bin/env python3
from __future__ import annotations

import gzip
import json
import math
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans


ROOT = Path("/Users/jijiabiao/Desktop/耳聋课题")
OUT = ROOT / "2026_supplement_analysis"
TEMP = Path("/tmp/ear_study_learning")
EXPRESSION = ROOT / "2022耳聋耳鸣" / "public" / "exp.txt"


def zscore_rows(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.sub(frame.mean(axis=1), axis=0).div(
        frame.std(axis=1).replace(0, np.nan),
        axis=0,
    )


def group_means(frame: pd.DataFrame, groups: dict[str, list[str]]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            group: frame[columns].mean(axis=1)
            for group, columns in groups.items()
        }
    )


def ordinal_correlation(
    frame: pd.DataFrame,
    values: np.ndarray,
    genes: list[str],
) -> tuple[float, float, float]:
    genes = [gene for gene in genes if gene in frame.index]
    if len(genes) < 3:
        return math.nan, math.nan, math.nan
    matrix = frame.loc[genes].to_numpy(dtype=float)
    matrix = (matrix - matrix.mean(axis=1, keepdims=True)) / np.where(
        matrix.std(axis=1, keepdims=True) == 0,
        1,
        matrix.std(axis=1, keepdims=True),
    )
    score = matrix.mean(axis=0)
    x = values - values.mean()
    y = score - score.mean()
    denominator = math.sqrt(float((x**2).sum()) * float((y**2).sum()))
    r = 0.0 if denominator == 0 else float((x * y).sum() / denominator)
    if len(values) <= 2 or abs(r) >= 1:
        p_value = 1.0 if abs(r) < 1 else 0.0
    else:
        t_value = r * math.sqrt((len(values) - 2) / max(1 - r**2, 1e-300))
        p_value = math.erfc(abs(t_value) / math.sqrt(2))
    return r, p_value, float(np.nanmean(matrix @ matrix.T))


def permutation_z(
    frame: pd.DataFrame,
    values: np.ndarray,
    genes: list[str],
    iterations: int = 1000,
    seed: int = 2026,
) -> tuple[float, float]:
    observed, _, _ = ordinal_correlation(frame, values, genes)
    if not math.isfinite(observed):
        return math.nan, math.nan
    all_genes = list(frame.index)
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(iterations):
        draw = rng.choice(all_genes, size=min(len(genes), len(all_genes)), replace=False)
        value, _, _ = ordinal_correlation(frame, values, list(draw))
        null.append(value)
    null_values = np.asarray(
        [value for value in null if math.isfinite(value)],
        dtype=float,
    )
    if len(null_values) < 20 or null_values.std(ddof=1) == 0:
        return math.nan, math.nan
    z_score = (observed - null_values.mean()) / max(null_values.std(ddof=1), 1e-9)
    empirical_p = (np.sum(np.abs(null_values) >= abs(observed)) + 1) / (iterations + 1)
    return float(z_score), float(empirical_p)


def load_gse49543() -> tuple[pd.DataFrame, dict[str, list[str]]]:
    expression = pd.read_csv(EXPRESSION, sep="\t", index_col=0)
    expression = expression.groupby(level=0).mean()
    groups = {
        group: [column for column in expression.columns if column.startswith(group + "-")]
        for group in ["YC", "MA", "MP", "SP"]
    }
    return np.log2(expression + 1.0), groups


def load_gse49522() -> tuple[pd.DataFrame, dict[str, list[str]]]:
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
        raise RuntimeError("Missing GSE49522 hearing status")
    labels = [value.strip().strip('"').split(":", 1)[-1].strip() for value in status_line.split("\t")[1:]]
    mapping = {
        "Young Control": "YC",
        "Middle-Aged": "MA",
        "Mild Presbycusis": "MP",
        "Severe Presbycusis": "SP",
    }
    groups = {"YC": [], "MA": [], "MP": [], "SP": []}
    for index, label in enumerate(labels):
        groups[mapping[label]].append(f"{mapping[label]}-{index + 1}")

    rows = []
    in_platform = False
    with (TEMP / "GPL339_full.txt").open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_platform = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_platform:
                rows.append(line.rstrip("\n"))
    platform = pd.read_csv(StringIO("\n".join(rows)), sep="\t", dtype=str)
    mapping_frame = platform[["ID", "Gene Symbol"]].dropna()
    mapping_frame = mapping_frame[
        (mapping_frame["Gene Symbol"] != "")
        & (mapping_frame["Gene Symbol"] != "---")
        & (~mapping_frame["Gene Symbol"].str.contains("///", regex=False))
    ]
    probe_to_symbol = dict(zip(mapping_frame["ID"], mapping_frame["Gene Symbol"]))
    matrix = pd.read_csv(StringIO("\n".join(lines)), sep="\t", index_col=0)
    ordered_names = [name for group in ("YC", "MA", "MP", "SP") for name in groups[group]]
    matrix.columns = ordered_names
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
    expression = np.log2(frame[old + young] + 1.0)
    age = np.array([1.0] * len(old) + [0.0] * len(young))
    return expression, age


def load_gse154833() -> tuple[pd.DataFrame, np.ndarray]:
    frame = pd.read_excel(
        TEMP / "GSE154833_Stria_Vascularis_RPKM.xlsx",
        sheet_name="Means",
    )
    frame = frame.dropna(subset=["Gene name"]).copy()
    frame["Gene name"] = frame["Gene name"].astype(str)
    columns = [
        "1mo_SV means",
        "9mo_SV means",
        "26mo_SV means",
    ]
    frame[columns] = frame[columns].apply(pd.to_numeric, errors="coerce")
    frame = frame.groupby("Gene name")[columns].mean()
    age = np.array([0.0, 1.0, 2.0])
    return frame, age


def main() -> int:
    expression, groups = load_gse49543()
    severity_genes = pd.read_csv(OUT / "severity_trend_all_genes.csv")
    severity_genes = (
        severity_genes[severity_genes["fdr_bh"] <= 0.05]
        .sort_values("fdr_bh")
        .head(200)
    )
    genes = [gene for gene in severity_genes["gene"].astype(str) if gene in expression.index]
    trajectory = group_means(expression, groups).loc[genes]
    trajectory = zscore_rows(trajectory).dropna()
    kmeans = KMeans(n_clusters=4, n_init=100, random_state=2026)
    cluster = kmeans.fit_predict(trajectory.to_numpy(dtype=float))
    assignment = pd.DataFrame(
        {
            "gene": trajectory.index,
            "module": [f"M{value + 1}" for value in cluster],
        }
    )
    assignment.to_csv(OUT / "module_assignment.csv", index=False, encoding="utf-8-sig")

    cochlea = expression
    central, central_groups = load_gse49522()
    external, external_age = load_gse233798()
    stria, stria_age = load_gse154833()

    cochlea_values = np.array([{"YC": 0, "MA": 1, "MP": 2, "SP": 3}[c.split("-")[0]] for c in cochlea.columns])
    central_values = np.array(
        [{"YC": 0, "MA": 1, "MP": 2, "SP": 3}[c.split("-")[0]] for c in central.columns]
    )

    rows = []
    for module, members in assignment.groupby("module")["gene"]:
        module_genes = [gene for gene in members.tolist() if gene in cochlea.index and gene in central.index]
        if len(module_genes) < 3:
            continue
        cohort_results = {}
        for dataset_name, frame, values in [
            ("cochlea_GSE49543", cochlea, cochlea_values),
            ("inferior_colliculus_GSE49522", central, central_values),
            ("cochlea_GSE233798", external, external_age),
            ("stria_vascularis_GSE154833", stria, stria_age),
        ]:
            available = [gene for gene in module_genes if gene in frame.index]
            r, p_value, coherence = ordinal_correlation(frame, values, available)
            z_score, empirical_p = permutation_z(frame, values, available)
            cohort_results[dataset_name] = {
                "genes": len(available),
                "correlation": r,
                "p_value": p_value,
                "coherence": coherence,
                "permutation_z": z_score,
                "permutation_p": empirical_p,
            }
        rows.append(
            {
                "module": module,
                "genes": len(module_genes),
                **{
                    f"{dataset}_{metric}": value
                    for dataset, metrics in cohort_results.items()
                    for metric, value in metrics.items()
                },
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "module_preservation_summary.csv", index=False, encoding="utf-8-sig")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
