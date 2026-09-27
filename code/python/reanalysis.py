#!/usr/bin/env python3
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path("/Users/jijiabiao/Desktop/耳聋课题")
OUT = ROOT / "2026_supplement_analysis"
TEMP = Path("/tmp/ear_study_learning")
EXPRESSION = ROOT / "2022耳聋耳鸣" / "public" / "exp.txt"

CONTRASTS = [
    ("MA_vs_YC", "MA", "YC"),
    ("MP_vs_MA", "MP", "MA"),
    ("MP_vs_YC", "MP", "YC"),
    ("SP_vs_MA", "SP", "MA"),
    ("SP_vs_MP", "SP", "MP"),
    ("SP_vs_YC", "SP", "YC"),
]

CANDIDATES = [
    "Lamp3",
    "Tbk1",
    "Ddx1",
    "Fgf16",
    "Setd1a",
    "Clec4d",
    "Clec7a",
    "Ctss",
    "Mpeg1",
    "Fcgr3",
    "Cd68",
    "Lgals3",
    "Laptm5",
    "C3ar1",
    "Ms4a7",
    "Cd74",
    "Mif",
    "Apoe",
    "C1qa",
    "C1qb",
    "Tyrobp",
    "Pik3cd",
    "Slc11a1",
    "Nfatc2",
    "H2-Aa",
    "H2-Eb1",
]


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    p_values = np.where(np.isfinite(p_values), p_values, 1.0)
    p_values = np.clip(p_values, np.nextafter(0.0, 1.0), 1.0)
    count = len(p_values)
    order = np.argsort(p_values)
    adjusted = p_values[order] * count / (np.arange(count) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    result = np.empty(count, dtype=float)
    result[order] = adjusted
    return result


def welch_ttest(first: np.ndarray, second: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean_first = first.mean(axis=1)
    mean_second = second.mean(axis=1)
    variance_first = first.var(axis=1, ddof=1)
    variance_second = second.var(axis=1, ddof=1)
    standard_error = np.sqrt(
        variance_first / first.shape[1] + variance_second / second.shape[1]
    )
    delta = mean_first - mean_second
    t_value = np.divide(
        delta,
        standard_error,
        out=np.zeros_like(delta, dtype=float),
        where=standard_error > 0,
    )
    p_value = np.array(
        [math.erfc(abs(value) / math.sqrt(2)) for value in t_value],
        dtype=float,
    )
    p_value = np.clip(p_value, np.nextafter(0.0, 1.0), 1.0)
    return t_value, p_value


def expression_pca(expression: pd.DataFrame) -> dict[str, Any]:
    sample_centered = expression.sub(expression.mean(axis=1), axis=0)
    matrix = sample_centered.to_numpy(dtype=float)
    _, singular_values, vt = np.linalg.svd(matrix, full_matrices=False)
    explained = singular_values**2
    explained = explained / explained.sum()
    scores = vt[:4].T * singular_values[:4]
    scores_frame = pd.DataFrame(
        scores,
        index=expression.columns,
        columns=["PC1", "PC2", "PC3", "PC4"],
    )
    scores_frame.insert(0, "group", [name.split("-")[0] for name in scores_frame.index])
    scores_frame.to_csv(OUT / "gse49543_pca_scores.csv", encoding="utf-8-sig")
    return {
        "explained_variance": [float(value) for value in explained[:10]],
    }


def severity_trend(expression: pd.DataFrame) -> tuple[dict[str, Any], pd.DataFrame]:
    severity = {"YC": 0.0, "MA": 1.0, "MP": 2.0, "SP": 3.0}
    values = np.array(
        [severity[column.split("-")[0]] for column in expression.columns],
        dtype=float,
    )
    matrix = expression.to_numpy(dtype=float)
    centered_severity = values - values.mean()
    centered_expression = matrix - matrix.mean(axis=1, keepdims=True)
    denominator = np.sqrt(
        np.square(centered_expression).sum(axis=1)
        * np.square(centered_severity).sum()
    )
    correlation = np.divide(
        centered_expression @ centered_severity,
        denominator,
        out=np.zeros(matrix.shape[0], dtype=float),
        where=denominator > 0,
    )
    correlation = np.clip(correlation, -1.0, 1.0)
    sample_count = len(values)
    t_value = correlation * np.sqrt(
        (sample_count - 2) / np.maximum(1.0 - np.square(correlation), 1e-300)
    )
    p_value = np.array(
        [math.erfc(abs(value) / math.sqrt(2)) for value in t_value],
        dtype=float,
    )
    p_value = np.clip(p_value, np.nextafter(0.0, 1.0), 1.0)
    fdr = benjamini_hochberg(p_value)
    frame = pd.DataFrame(
        {
            "gene": expression.index,
            "severity_correlation": correlation,
            "welch_t": t_value,
            "p_value": p_value,
            "fdr_bh": fdr,
        }
    ).set_index("gene")
    frame.sort_values(["fdr_bh", "severity_correlation"], ascending=[True, False]).to_csv(
        OUT / "severity_trend_all_genes.csv",
        encoding="utf-8-sig",
    )
    significant = frame[frame["fdr_bh"] <= 0.05]
    summary = {
        "significant_genes": len(significant),
        "positive_genes": int((significant["severity_correlation"] > 0).sum()),
        "negative_genes": int((significant["severity_correlation"] < 0).sum()),
        "top_positive": significant[significant["severity_correlation"] > 0]
        .sort_values(["fdr_bh", "severity_correlation"], ascending=[True, False])
        .head(50)
        .reset_index()
        .to_dict(orient="records"),
        "top_negative": significant[significant["severity_correlation"] < 0]
        .sort_values(["fdr_bh", "severity_correlation"], ascending=[True, True])
        .head(50)
        .reset_index()
        .to_dict(orient="records"),
    }
    return summary, frame


def analyse_gse49543() -> tuple[dict[str, Any], dict[str, pd.DataFrame]]:
    expression = pd.read_csv(EXPRESSION, sep="\t", index_col=0)
    expression.index = expression.index.astype(str)
    expression = expression.groupby(level=0).mean()
    logged = np.log2(expression + 1.0)

    result_tables: dict[str, pd.DataFrame] = {}
    summary: dict[str, Any] = {}
    correlation = logged.corr(method="pearson")
    groups = {
        group: [column for column in logged.columns if column.startswith(group + "-")]
        for group in ["YC", "MA", "MP", "SP"]
    }
    within_group = {}
    for group, columns in groups.items():
        block = correlation.loc[columns, columns].to_numpy(dtype=float)
        values = block[np.triu_indices(len(columns), 1)]
        within_group[group] = float(np.nanmean(values))
    summary["samples"] = {
        group: columns for group, columns in groups.items()
    }
    summary["within_group_correlation"] = within_group
    summary["pca"] = expression_pca(logged)

    for name, first_group, second_group in CONTRASTS:
        first = logged[groups[first_group]].to_numpy(dtype=float)
        second = logged[groups[second_group]].to_numpy(dtype=float)
        t_value, p_value = welch_ttest(first, second)
        fdr = benjamini_hochberg(p_value)
        log2_fold_change = first.mean(axis=1) - second.mean(axis=1)
        table = pd.DataFrame(
            {
                "gene": logged.index,
                "mean_first_log2": first.mean(axis=1),
                "mean_second_log2": second.mean(axis=1),
                "log2FC": log2_fold_change,
                "fold_change": np.power(2.0, log2_fold_change),
                "welch_t": t_value,
                "p_value": p_value,
                "fdr_bh": fdr,
            }
        ).set_index("gene")
        table = table.sort_values(["fdr_bh", "p_value", "log2FC"], ascending=[True, True, False])
        table.to_csv(
            OUT / f"deg_{name}_log2_welch_bh.csv",
            encoding="utf-8-sig",
        )
        result_tables[name] = table
        summary[name] = {
            "counts": {
                "fdr_0.05_log2fc_0.585": {
                    "up": int(((table["fdr_bh"] <= 0.05) & (table["log2FC"] >= 0.585)).sum()),
                    "down": int(((table["fdr_bh"] <= 0.05) & (table["log2FC"] <= -0.585)).sum()),
                },
                "fdr_0.05_log2fc_1": {
                    "up": int(((table["fdr_bh"] <= 0.05) & (table["log2FC"] >= 1.0)).sum()),
                    "down": int(((table["fdr_bh"] <= 0.05) & (table["log2FC"] <= -1.0)).sum()),
                },
            },
            "top_up": table[table["log2FC"] > 0]
            .sort_values(["fdr_bh", "log2FC"], ascending=[True, False])
            .head(20)
            .reset_index()
            .to_dict(orient="records"),
            "top_down": table[table["log2FC"] < 0]
            .sort_values(["fdr_bh", "log2FC"], ascending=[True, True])
            .head(20)
            .reset_index()
            .to_dict(orient="records"),
        }

    rows = []
    for name, first_group, second_group in CONTRASTS:
        table = result_tables[name]
        for gene in CANDIDATES:
            if gene not in table.index:
                continue
            row = table.loc[gene]
            rows.append(
                {
                    "gene": gene,
                    "contrast": name,
                    "first_group": first_group,
                    "second_group": second_group,
                    "log2FC": float(row["log2FC"]),
                    "fold_change": float(row["fold_change"]),
                    "p_value": float(row["p_value"]),
                    "fdr_bh": float(row["fdr_bh"]),
                }
            )
    candidate_frame = pd.DataFrame(rows)
    candidate_frame.to_csv(
        OUT / "candidate_gene_effects_gse49543.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary["candidate_genes"] = rows
    return summary, result_tables


def read_gse233798() -> tuple[pd.DataFrame, list[str], list[str]]:
    path = TEMP / "GSE233798_FPKM.xlsx"
    frame = pd.read_excel(path)
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
    return frame, young, old


def standardise_validation_frame(
    frame: pd.DataFrame,
    group_one: list[str],
    group_two: list[str],
    dataset: str,
) -> pd.DataFrame:
    first = np.log2(frame[group_one].to_numpy(dtype=float) + 1.0)
    second = np.log2(frame[group_two].to_numpy(dtype=float) + 1.0)
    t_value, p_value = welch_ttest(first, second)
    delta = first.mean(axis=1) - second.mean(axis=1)
    fdr = benjamini_hochberg(p_value)
    return pd.DataFrame(
        {
            "gene": frame.index,
            "dataset": dataset,
            "log2FC_target_minus_control": delta,
            "fold_change": np.power(2.0, delta),
            "welch_t": t_value,
            "p_value": p_value,
            "fdr_bh": fdr,
        }
    ).set_index("gene")


def validate_gse233798() -> tuple[pd.DataFrame, dict[str, Any]]:
    frame, young, old = read_gse233798()
    result = standardise_validation_frame(frame, old, young, "GSE233798_old_vs_young")
    result.to_csv(OUT / "validation_GSE233798_all_genes.csv", encoding="utf-8-sig")
    candidates = result.reindex(CANDIDATES).dropna(subset=["log2FC_target_minus_control"])
    candidates.reset_index().to_csv(
        OUT / "validation_GSE233798_candidate_genes.csv",
        index=False,
        encoding="utf-8-sig",
    )
    return result, {
        "samples": {"young": young, "old": old},
        "genes": len(result),
        "candidate_genes": candidates.reset_index().to_dict(orient="records"),
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    gse49543_summary, _ = analyse_gse49543()
    expression = pd.read_csv(EXPRESSION, sep="\t", index_col=0)
    expression = expression.groupby(level=0).mean()
    logged = np.log2(expression + 1.0)
    severity_summary, _ = severity_trend(logged)
    gse49543_summary["severity_trend"] = severity_summary
    validation_summary: dict[str, Any] = {}
    try:
        _, validation_summary["GSE233798"] = validate_gse233798()
    except Exception as exc:
        validation_summary["GSE233798"] = {"error": str(exc)}

    summary = {
        "dataset": "GSE49543",
        "notes": [
            "GSE49543 is an Affymetrix MOE430A microarray dataset, not RNA-seq.",
            "Values were transformed with log2(x + 1), then Welch t tests and BH FDR correction were applied.",
            "This is a conservative sensitivity analysis and is not a replacement for a limma/eBayes workflow.",
        ],
        "gse49543": gse49543_summary,
        "validation": validation_summary,
    }
    import json

    (OUT / "supplement_analysis_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )
    print(json.dumps(
        {
            "output_dir": str(OUT),
            "gse49543_counts": {
                key: value["counts"]
                for key, value in gse49543_summary.items()
                if isinstance(value, dict) and "counts" in value
            },
            "validation": {
                key: value.get("genes", "error")
                if isinstance(value, dict)
                else None
                for key, value in validation_summary.items()
            },
        },
        ensure_ascii=False,
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
