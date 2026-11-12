#!/usr/bin/env python3
"""Prepare the three external GEO expression matrices for Welch/BH analysis.

The GEO submitters supplied normalized FPKM or RPKM values. This script keeps
that normalization, applies log2(x + 1), resolves duplicate gene symbols, and
records the sample and run mapping without treating GSE154833 lanes as
independent biological replicates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import openpyxl


DATASETS = {
    "GSE233798": {
        "file": "GSE233798_Processed_data_files_FPKM_.xlsx",
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/"
            "GSE233798/suppl/GSE233798_Processed_data_files_FPKM_.xlsx"
        ),
        "unit": "FPKM",
    },
    "GSE153882": {
        "file": "GSE153882_Aging_Hair_Cells_RPKM.xlsx",
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE153nnn/"
            "GSE153882/suppl/GSE153882_Aging_Hair_Cells_RPKM.xlsx"
        ),
        "unit": "RPKM",
    },
    "GSE154833": {
        "file": "GSE154833_Stria_Vascularis_RPKM.xlsx",
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/"
            "GSE154833/suppl/GSE154833_Stria_Vascularis_RPKM.xlsx"
        ),
        "unit": "RPKM",
    },
}

GSE233798_MANIFEST = [
    ("GSE233798", "cochlea", "whole_cochlea", 0, "young", "young1", "GSM7437033", 1),
    ("GSE233798", "cochlea", "whole_cochlea", 0, "young", "young2", "GSM7437034", 1),
    ("GSE233798", "cochlea", "whole_cochlea", 0, "young", "young3", "GSM7437035", 1),
    ("GSE233798", "cochlea", "whole_cochlea", 1, "old", "old1", "GSM7437036", 1),
    ("GSE233798", "cochlea", "whole_cochlea", 1, "old", "old2", "GSM7437037", 1),
    ("GSE233798", "cochlea", "whole_cochlea", 1, "old", "old3", "GSM7437038", 1),
]

GSE153882_MANIFEST = [
    ("GSE153882", "cochlea", "IHC", 9, "9m", "9m_IHC_1", "GSM5595453", 1),
    ("GSE153882", "cochlea", "IHC", 9, "9m", "9m_IHC_2", "GSM5595454", 1),
    ("GSE153882", "cochlea", "IHC", 9, "9m", "9m_IHC_3", "GSM5595455", 1),
    ("GSE153882", "cochlea", "IHC", 9, "9m", "9m_IHC_4", "GSM5595456", 1),
    ("GSE153882", "cochlea", "OHC", 9, "9m", "9m_OHC_1", "GSM5595457", 1),
    ("GSE153882", "cochlea", "OHC", 9, "9m", "9m_OHC_2", "GSM5595458", 1),
    ("GSE153882", "cochlea", "OHC", 9, "9m", "9m_OHC_3", "GSM5595459", 1),
    ("GSE153882", "cochlea", "OHC", 9, "9m", "9m_OHC_4", "GSM5595460", 1),
    ("GSE153882", "cochlea", "IHC", 26, "26m", "26m_IHC_1", "GSM5595461", 1),
    ("GSE153882", "cochlea", "IHC", 26, "26m", "26m_IHC_2", "GSM5595462", 1),
    ("GSE153882", "cochlea", "IHC", 26, "26m", "26m_IHC_3", "GSM5595463", 1),
    ("GSE153882", "cochlea", "IHC", 26, "26m", "26m_IHC_4", "GSM5595464", 1),
    ("GSE153882", "cochlea", "OHC", 26, "26m", "26m_OHC_1", "GSM5595465", 1),
    ("GSE153882", "cochlea", "OHC", 26, "26m", "26m_OHC_2", "GSM5595466", 1),
    ("GSE153882", "cochlea", "OHC", 26, "26m", "26m_OHC_3", "GSM5595467", 1),
    ("GSE153882", "cochlea", "OHC", 26, "26m", "26m_OHC_4", "GSM5595468", 1),
]

GSE154833_RUNS = [
    ("GSM4681358", "CBA_P30_SV_A_L001"),
    ("GSM4681359", "CBA_P30_SV_A_L002"),
    ("GSM4681360", "CBA_P30_SV_A_L003"),
    ("GSM4681361", "CBA_P30_SV_A_L004"),
    ("GSM4681362", "CBA_P30_SV_B_L001"),
    ("GSM4681363", "CBA_P30_SV_B_L002"),
    ("GSM4681364", "CBA_P30_SV_B_L003"),
    ("GSM4681365", "CBA_P30_SV_B_L004"),
    ("GSM4681366", "CBA_P30_SV_C_L001"),
    ("GSM4681367", "CBA_P30_SV_C_L002"),
    ("GSM4681368", "CBA_P30_SV_C_L003"),
    ("GSM4681369", "CBA_P30_SV_C_L004"),
    ("GSM4681370", "CBA_9M_SV_A_L001"),
    ("GSM4681371", "CBA_9M_SV_A_L002"),
    ("GSM4681372", "CBA_9M_SV_A_L003"),
    ("GSM4681373", "CBA_9M_SV_A_L004"),
    ("GSM4681374", "CBA_9M_SV_B_L001"),
    ("GSM4681375", "CBA_9M_SV_B_L002"),
    ("GSM4681376", "CBA_9M_SV_B_L003"),
    ("GSM4681377", "CBA_9M_SV_B_L004"),
    ("GSM4681378", "CBA_9M_SV_C_L001"),
    ("GSM4681379", "CBA_9M_SV_C_L002"),
    ("GSM4681380", "CBA_9M_SV_C_L003"),
    ("GSM4681381", "CBA_9M_SV_C_L004"),
    ("GSM4681382", "CBA_26M_SV_A_L001"),
    ("GSM4681383", "CBA_26M_SV_A_L002"),
    ("GSM4681384", "CBA_26M_SV_A_L003"),
    ("GSM4681385", "CBA_26M_SV_A_L004"),
    ("GSM4681386", "CBA_26M_SV_B_L001"),
    ("GSM4681387", "CBA_26M_SV_B_L002"),
    ("GSM4681388", "CBA_26M_SV_B_L003"),
    ("GSM4681389", "CBA_26M_SV_B_L004"),
    ("GSM4681390", "CBA_26M_SV_C_L001"),
    ("GSM4681391", "CBA_26M_SV_C_L002"),
    ("GSM4681392", "CBA_26M_SV_C_L003"),
    ("GSM4681393", "CBA_26M_SV_C_L004"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean_symbol(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    symbol = str(value).strip()
    if symbol.lower() in {"", "na", "nan", "none"}:
        return ""
    return symbol


def first_data_sheet(path: Path) -> str:
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    for name in workbook.sheetnames:
        sheet = workbook[name]
        for row in sheet.iter_rows(values_only=True):
            if any(value is not None for value in row):
                workbook.close()
                return name
    workbook.close()
    raise ValueError(f"No non-empty sheet found in {path}")


def read_excel_frame(path: Path, sheet_name: str | None = None) -> pd.DataFrame:
    sheet = sheet_name or first_data_sheet(path)
    frame = pd.read_excel(path, sheet_name=sheet, engine="openpyxl")
    frame.columns = [str(column).strip() for column in frame.columns]
    return frame


def collapse_symbols(
    frame: pd.DataFrame,
    symbol_column: str,
    sample_columns: list[str],
    dataset: str,
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    if symbol_column not in frame.columns:
        raise ValueError(f"{dataset}: missing symbol column {symbol_column!r}")
    frame = frame.copy()
    source_rows = len(frame)
    frame["gene"] = frame[symbol_column].map(clean_symbol)
    for column in sample_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    missing_value_count = int(frame[sample_columns].isna().sum().sum())
    removed_missing_symbol = int((frame["gene"] == "").sum())
    frame = frame.loc[frame["gene"] != ""].copy()
    genes_after_symbol_filter = len(frame)
    incomplete_rows = frame[sample_columns].isna().any(axis=1)
    rows_excluded_incomplete = int(incomplete_rows.sum())
    frame = frame.loc[~incomplete_rows].copy()
    frame["_source_mean"] = frame[sample_columns].mean(axis=1, skipna=True)
    duplicated = frame["gene"].duplicated(keep=False)
    duplicate_rows = frame.loc[duplicated, ["gene", symbol_column, "_source_mean"]].copy()
    duplicate_rows = duplicate_rows.sort_values(
        ["gene", "_source_mean"], ascending=[True, False]
    )
    duplicate_symbols = int(duplicate_rows["gene"].nunique())
    duplicate_rows_discarded = int(duplicated.sum() - duplicate_symbols)
    frame = frame.sort_values(
        ["gene", "_source_mean"], ascending=[True, False]
    ).drop_duplicates("gene", keep="first")
    frame = frame.loc[:, ["gene", *sample_columns]].copy()
    frame = frame.set_index("gene", drop=True)
    collapse_log = [
        {
            "dataset": dataset,
            "source_rows": source_rows,
            "genes_after_symbol_filter": genes_after_symbol_filter,
            "rows_removed_missing_symbol": removed_missing_symbol,
            "rows_excluded_incomplete_sample_values": rows_excluded_incomplete,
            "duplicate_symbol_rows_discarded": duplicate_rows_discarded,
            "genes_with_duplicate_symbol_rows": duplicate_symbols,
            "missing_sample_values": missing_value_count,
            "gene_exclusion_rule": (
                "exclude rows with a missing or non-numeric value in any selected "
                "biological-sample column"
            ),
            "collapse_rule": "keep symbol row with highest mean submitted normalized value",
        }
    ]
    return frame, collapse_log


def matrix_to_frames(matrix: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = matrix.reset_index()
    logged = np.log2(matrix.astype(float) + 1.0)
    log_frame = logged.reset_index()
    return raw, log_frame


def prepare_gse233798(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    frame = read_excel_frame(path, "gene_expression（1）")
    sample_columns = [
        "fpkm_Base_young1",
        "fpkm_Base_young2",
        "fpkm_Base_young3",
        "fpkm_Base_old1",
        "fpkm_Base_old2",
        "fpkm_Base_old3",
    ]
    matrix, collapse_log = collapse_symbols(
        frame, "gene_symbol", sample_columns, "GSE233798"
    )
    return (*matrix_to_frames(matrix), collapse_log)


def prepare_gse153882(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    sample_pattern = re.compile(
        r"(?:9m|26m)_(?:IHCs|OHCs)_[1-4](?:\s+RPKM value)?",
        flags=re.IGNORECASE,
    )
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    selected_sheet = None
    for sheet_name in workbook.sheetnames:
        sheet = workbook[sheet_name]
        rows = sheet.iter_rows(values_only=True)
        for row in rows:
            cleaned = [str(value).strip() if value is not None else "" for value in row]
            if sum(bool(value) for value in cleaned) < 4:
                continue
            matches = sum(
                bool(sample_pattern.fullmatch(value))
                for value in cleaned
            )
            if matches >= 8:
                selected_sheet = sheet_name
                break
        if selected_sheet:
            break
    workbook.close()
    if not selected_sheet:
        raise ValueError("GSE153882: could not locate a sheet with 16 sample columns")
    frame = read_excel_frame(path, selected_sheet)
    sample_columns = [
        column
        for column in frame.columns
        if sample_pattern.fullmatch(column)
    ]
    if len(sample_columns) != 16:
        raise ValueError(
            f"GSE153882: expected 16 sample columns, found {len(sample_columns)}"
        )
    rename_map = {
        column: re.match(r"((?:9m|26m)_(?:IHCs|OHCs)_[1-4])", column).group(1)
        .replace("IHCs", "IHC")
        .replace("OHCs", "OHC")
        for column in sample_columns
    }
    frame = frame.rename(columns=rename_map)
    sample_columns = [rename_map[column] for column in sample_columns]
    symbol_column = next(
        (
            column
            for column in frame.columns
            if column.lower() in {"gene_name", "gene symbol", "gene_symbol", "gene"}
        ),
        None,
    )
    if symbol_column is None:
        symbol_column = frame.columns[0]
    matrix, collapse_log = collapse_symbols(
        frame, symbol_column, sample_columns, "GSE153882"
    )
    return (*matrix_to_frames(matrix), collapse_log)


def prepare_gse154833(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    frame = read_excel_frame(path, "All values")
    sample_columns = [
        column
        for column in frame.columns
        if re.fullmatch(r"(?:1mo|9mo|26mo)_SV_[ABC]", column)
    ]
    expected = [f"{age}_SV_{rep}" for age in ("1mo", "9mo", "26mo") for rep in "ABC"]
    if sorted(sample_columns) != sorted(expected):
        raise ValueError(
            "GSE154833: expected nine biological-sample columns "
            f"{expected}; found {sample_columns}"
        )
    matrix, collapse_log = collapse_symbols(
        frame, "Gene name", sample_columns, "GSE154833"
    )
    return (*matrix_to_frames(matrix), collapse_log)


def build_biological_sample_manifest() -> pd.DataFrame:
    rows = []
    for dataset, tissue, cell_type, age, group, sample_id, gsm, runs in (
        GSE233798_MANIFEST + GSE153882_MANIFEST
    ):
        rows.append(
            {
                "dataset": dataset,
                "tissue": tissue,
                "cell_type": cell_type,
                "age_months": age,
                "group": group,
                "sample_id": sample_id,
                "gsm_ids": gsm,
                "n_gsm_runs": runs,
                "replicate": sample_id.rsplit("_", 1)[-1],
                "analysis_level": "biological_sample",
            }
        )
    for age, age_token, group in (
        (1, "1mo", "1m"),
        (9, "9mo", "9m"),
        (26, "26mo", "26m"),
    ):
        geo_age_token = {1: "P30", 9: "9M", 26: "26M"}[age]
        for replicate in "ABC":
            sample_id = f"{age_token}_SV_{replicate}"
            gsm_ids = [
                gsm
                for gsm, title in GSE154833_RUNS
                if re.fullmatch(
                    rf"CBA_{geo_age_token}_SV_{replicate}_L00[1-4]",
                    title,
                )
            ]
            rows.append(
                {
                    "dataset": "GSE154833",
                    "tissue": "stria vascularis",
                    "cell_type": "stria vascularis",
                    "age_months": age,
                    "group": group,
                    "sample_id": sample_id,
                    "gsm_ids": ";".join(gsm_ids),
                    "n_gsm_runs": len(gsm_ids),
                    "replicate": replicate,
                    "analysis_level": "biological_sample",
                }
            )
    return pd.DataFrame(rows)


def build_run_manifest() -> pd.DataFrame:
    rows = []
    for dataset, _, cell_type, age, group, sample_id, gsm, _ in (
        GSE233798_MANIFEST + GSE153882_MANIFEST
    ):
        rows.append(
            {
                "dataset": dataset,
                "gsm_id": gsm,
                "gsm_title": sample_id,
                "biological_sample_id": sample_id,
                "cell_type": cell_type,
                "age_months": age,
                "group": group,
                "lane": "",
                "included": "yes",
            }
        )
    for gsm, title in GSE154833_RUNS:
        match = re.fullmatch(
            r"CBA_(P30|9M|26M)_SV_([ABC])_(L00[1-4])", title
        )
        if not match:
            raise ValueError(f"Unexpected GSE154833 sample title: {title}")
        age_token, replicate, lane = match.groups()
        age = {"P30": 1, "9M": 9, "26M": 26}[age_token]
        group = {1: "1m", 9: "9m", 26: "26m"}[age]
        workbook_token = {1: "1mo", 9: "9mo", 26: "26mo"}[age]
        rows.append(
            {
                "dataset": "GSE154833",
                "gsm_id": gsm,
                "gsm_title": title,
                "biological_sample_id": f"{workbook_token}_SV_{replicate}",
                "cell_type": "stria vascularis",
                "age_months": age,
                "group": group,
                "lane": lane,
                "included": "yes",
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--archive-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    args = parser.parse_args()
    root = args.archive_root.resolve()
    input_dir = root / "inputs" / "external" / "geo"
    output_dir = root / "inputs" / "external" / "prepared"
    output_dir.mkdir(parents=True, exist_ok=True)

    prepared: dict[str, tuple[pd.DataFrame, pd.DataFrame, list[dict]]] = {}
    for accession in ("GSE233798", "GSE153882", "GSE154833"):
        path = input_dir / DATASETS[accession]["file"]
        if not path.exists():
            raise FileNotFoundError(path)
        if accession == "GSE233798":
            prepared[accession] = prepare_gse233798(path)
        elif accession == "GSE153882":
            prepared[accession] = prepare_gse153882(path)
        else:
            prepared[accession] = prepare_gse154833(path)

        raw, logged, _ = prepared[accession]
        raw.to_csv(output_dir / f"{accession}_expression_raw.csv", index=False)
        logged.to_csv(output_dir / f"{accession}_expression_log2.csv", index=False)

    samples = build_biological_sample_manifest()
    runs = build_run_manifest()
    samples.to_csv(output_dir / "sample_mapping_external.csv", index=False)
    runs.to_csv(output_dir / "run_mapping_external.csv", index=False)

    collapse_log = pd.DataFrame(
        [row for _, _, rows in prepared.values() for row in rows]
    )
    collapse_log.to_csv(output_dir / "gene_symbol_collapse_log.csv", index=False)

    exclusion_rows = [
            {
                "dataset": accession,
                "scope": "sample",
                "item_id": "all_mapped_samples",
                "decision": "retained",
                "n_items": int((samples["dataset"] == accession).sum()),
                "reason": "no additional sample exclusion was prespecified or applied",
                "source": "GEO series metadata and submitted processed workbook",
            }
            for accession in ("GSE233798", "GSE153882", "GSE154833")
        ]
    for row in collapse_log.to_dict(orient="records"):
        exclusion_rows.append(
            {
                "dataset": row["dataset"],
                "scope": "gene",
                "item_id": "rows_with_missing_sample_values",
                "decision": "excluded",
                "n_items": row["rows_excluded_incomplete_sample_values"],
                "reason": row["gene_exclusion_rule"],
                "source": "submitted processed expression workbook",
            }
        )
        exclusion_rows.append(
            {
                "dataset": row["dataset"],
                "scope": "gene_symbol",
                "item_id": "duplicate_symbol_rows",
                "decision": "collapsed",
                "n_items": row["duplicate_symbol_rows_discarded"],
                "reason": row["collapse_rule"],
                "source": "submitted processed expression workbook",
            }
        )
    exclusion_log = pd.DataFrame(exclusion_rows)
    exclusion_log.to_csv(output_dir / "exclusion_log_external.csv", index=False)

    external_results_dir = root / "results" / "external"
    output_files = sorted(
        [
            path
            for path in (
                list(output_dir.glob("*.csv"))
                + list(external_results_dir.glob("*"))
            )
            if path.is_file()
            and path.name != "external_input_manifest.json"
        ]
    )
    input_files = [input_dir / DATASETS[name]["file"] for name in DATASETS]
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script": "code/python/prepare_external_inputs.py",
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "openpyxl": openpyxl.__version__,
        "normalization": {
            "GSE233798": "submitted normalized FPKM; log2(FPKM + 1)",
            "GSE153882": "submitted normalized RPKM; log2(RPKM + 1)",
            "GSE154833": "submitted normalized RPKM; log2(RPKM + 1)",
        },
        "additional_cross_sample_normalization": "none",
        "duplicate_symbol_rule": (
            "retain the row with the highest mean submitted normalized value"
        ),
        "gse154833_design": (
            "36 GSM runs map to nine biological samples; the submitted All values "
            "sheet provides nine normalized biological-sample columns (A/B/C per age), "
            "so the analysis does not treat lanes as independent replicates"
        ),
        "sample_exclusions": "none beyond rows without a usable gene symbol",
        "inputs": [
            {
                "dataset": name,
                "path": str(path.relative_to(root)),
                "url": DATASETS[name]["url"],
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for name, path in zip(DATASETS, input_files)
        ],
        "outputs": [
            {
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in output_files
            if path.exists()
        ],
        "matrix_dimensions": {
            name: {
                "genes": int(len(logged)),
                "biological_samples": int(logged.shape[1] - 1),
            }
            for name, (_, logged, _) in prepared.items()
        },
        "duplicate_symbol_collapse": collapse_log.to_dict(orient="records"),
    }
    manifest_path = output_dir / "external_input_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest["matrix_dimensions"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
