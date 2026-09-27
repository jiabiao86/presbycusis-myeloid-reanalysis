#!/usr/bin/env python3
from __future__ import annotations

import gzip
import math
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/Users/jijiabiao/Desktop/耳聋课题")
OUT = ROOT / "2026_supplement_analysis"
TEMP = Path("/tmp/ear_study_learning")
MATRIX_GZ = TEMP / "GSE49522_series_matrix.txt.gz"
PLATFORM = TEMP / "GPL339_full.txt"


def bh(p_values: np.ndarray) -> np.ndarray:
    p_values = np.clip(np.asarray(p_values, dtype=float), np.nextafter(0.0, 1.0), 1.0)
    count = len(p_values)
    order = np.argsort(p_values)
    adjusted = p_values[order] * count / (np.arange(count) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(count, dtype=float)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def welch(first: np.ndarray, second: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    first_mean = first.mean(axis=1)
    second_mean = second.mean(axis=1)
    standard_error = np.sqrt(
        first.var(axis=1, ddof=1) / first.shape[1]
        + second.var(axis=1, ddof=1) / second.shape[1]
    )
    delta = first_mean - second_mean
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
    return t_value, np.clip(p_value, np.nextafter(0.0, 1.0), 1.0)


def severity_trend(expression: pd.DataFrame) -> pd.DataFrame:
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
    return pd.DataFrame(
        {
            "gene": expression.index,
            "central_severity_correlation": correlation,
            "central_p_value": p_value,
            "central_fdr_bh": bh(p_value),
        }
    ).set_index("gene")


def read_platform() -> pd.DataFrame:
    rows = []
    in_table = False
    with PLATFORM.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                rows.append(line.rstrip("\n"))
    frame = pd.read_csv(StringIO("\n".join(rows)), sep="\t", dtype=str)
    return frame


def read_matrix() -> tuple[pd.DataFrame, dict[str, list[str]]]:
    lines = []
    status_line = None
    in_table = False
    with gzip.open(MATRIX_GZ, "rt", encoding="utf-8", errors="replace") as handle:
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
        raise RuntimeError("Could not find hearing status metadata")
    labels = [
        value.strip().strip('"')
        for value in status_line.split("\t")[1:]
    ]
    grouping = {
        "Young Control": "YC",
        "Middle-Aged": "MA",
        "Mild Presbycusis": "MP",
        "Severe Presbycusis": "SP",
    }
    groups: dict[str, list[str]] = {"YC": [], "MA": [], "MP": [], "SP": []}
    for index, label in enumerate(labels):
        status = label.split(":", 1)[-1].strip()
        groups[grouping[status]].append(f"{grouping[status]}-{index + 1}")
    matrix = pd.read_csv(StringIO("\n".join(lines)), sep="\t", index_col=0)
    ordered_names = [name for group in ("YC", "MA", "MP", "SP") for name in groups[group]]
    if len(ordered_names) != len(matrix.columns):
        raise RuntimeError("Sample metadata and expression matrix columns do not match")
    matrix.columns = ordered_names
    return matrix, groups


def collapse_to_symbols(
    matrix: pd.DataFrame,
    platform: pd.DataFrame,
) -> pd.DataFrame:
    mapping = platform[["ID", "Gene Symbol"]].dropna()
    mapping = mapping[
        (mapping["Gene Symbol"] != "")
        & (mapping["Gene Symbol"] != "---")
        & (~mapping["Gene Symbol"].str.contains("///", regex=False))
    ]
    probe_to_symbol = dict(zip(mapping["ID"], mapping["Gene Symbol"]))
    symbol_frame = matrix.copy()
    symbol_frame.index = [probe_to_symbol.get(str(index), "") for index in symbol_frame.index]
    symbol_frame = symbol_frame[symbol_frame.index != ""]
    return symbol_frame.groupby(level=0).mean()


def main() -> int:
    matrix, groups = read_matrix()
    platform = read_platform()
    expression = collapse_to_symbols(matrix, platform)
    logged = np.log2(expression + 1.0)

    trend = severity_trend(logged)
    trend.sort_values("central_fdr_bh").to_csv(
        OUT / "gse49522_severity_trend.csv",
        encoding="utf-8-sig",
    )

    contrasts = [
        ("MA_vs_YC", "MA", "YC"),
        ("MP_vs_MA", "MP", "MA"),
        ("SP_vs_MA", "SP", "MA"),
        ("SP_vs_YC", "SP", "YC"),
    ]
    for name, first_group, second_group in contrasts:
        first = logged[groups[first_group]].to_numpy(dtype=float)
        second = logged[groups[second_group]].to_numpy(dtype=float)
        t_value, p_value = welch(first, second)
        result = pd.DataFrame(
            {
                "gene": logged.index,
                "mean_first_log2": first.mean(axis=1),
                "mean_second_log2": second.mean(axis=1),
                "log2FC": first.mean(axis=1) - second.mean(axis=1),
                "welch_t": t_value,
                "p_value": p_value,
                "fdr_bh": bh(p_value),
            }
        ).set_index("gene")
        result.sort_values("fdr_bh").to_csv(
            OUT / f"gse49522_deg_{name}.csv",
            encoding="utf-8-sig",
        )

    cochlea = pd.read_csv(OUT / "severity_trend_all_genes.csv").set_index("gene")
    common = cochlea.join(trend, how="inner")
    common["same_direction"] = (
        common["severity_correlation"] * common["central_severity_correlation"] > 0
    )
    common["shared_significant"] = (
        (common["fdr_bh"] <= 0.05) & (common["central_fdr_bh"] <= 0.05)
    )
    common.sort_values(["shared_significant", "central_fdr_bh"], ascending=[False, True]).to_csv(
        OUT / "peripheral_central_severity_comparison.csv",
        encoding="utf-8-sig",
    )

    candidates = [
        "Cd74",
        "H2-Aa",
        "H2-Eb1",
        "Setd1a",
        "Nfatc2",
        "Ctss",
        "Fcgr3",
        "Cd68",
        "Tyrobp",
        "Clec4d",
        "Mpeg1",
        "Lgals3",
        "C1qa",
        "C1qb",
        "C3ar1",
        "Ms4a7",
        "Apoe",
        "Mif",
        "Ddx1",
        "Lamp3",
        "Tbk1",
        "Fgf16",
    ]
    candidate = common.reindex(candidates).dropna(subset=["severity_correlation"]).reset_index()
    candidate.to_csv(
        OUT / "peripheral_central_candidate_comparison.csv",
        index=False,
        encoding="utf-8-sig",
    )

    x = common["severity_correlation"].to_numpy(dtype=float)
    y = common["central_severity_correlation"].to_numpy(dtype=float)
    correlation = float(np.corrcoef(x, y)[0, 1])
    if abs(correlation) >= 1.0:
        p_value = 0.0
    else:
        t_value = correlation * math.sqrt(
            (len(x) - 2) / max(1.0 - correlation**2, 1e-300)
        )
        p_value = math.erfc(abs(t_value) / math.sqrt(2))
    summary = {
        "central_genes": len(trend),
        "common_genes": len(common),
        "shared_significant_genes": int(common["shared_significant"].sum()),
        "same_direction_genes": int(common["same_direction"].sum()),
        "severity_correlation_pearson_r": float(correlation),
        "severity_correlation_p_value": float(p_value),
    }
    pd.Series(summary).to_json(
        OUT / "peripheral_central_summary.json",
        indent=2,
        force_ascii=False,
    )
    print(summary)
    print(candidate.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
