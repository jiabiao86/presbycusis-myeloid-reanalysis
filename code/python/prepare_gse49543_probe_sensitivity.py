#!/usr/bin/env python3
"""Prepare GSE49543 matrices under three probe-to-symbol aggregation rules."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_GEO_DIR = REPO_ROOT / "inputs" / "geo"
DEFAULT_OUTPUT = REPO_ROOT / "inputs" / "prepared" / "probe_sensitivity"
MONTHS = (
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def find_soffice() -> str:
    configured = os.environ.get("V14_SOFFICE", "").strip()
    if configured:
        path = shutil.which(configured) or configured
        if Path(path).is_file():
            return str(Path(path).resolve())
        raise FileNotFoundError(f"V14_SOFFICE does not resolve to a file: {configured}")
    located = shutil.which("soffice")
    if located:
        return located
    raise FileNotFoundError(
        "LibreOffice is required to convert the GSE49543 .xls workbook."
    )


def convert_workbook(source_gz: Path, work_dir: Path) -> Path:
    source_xls = work_dir / "GSE49543_further_processed.xls"
    with gzip.open(source_gz, "rb") as source, source_xls.open("wb") as target:
        shutil.copyfileobj(source, target)
    profile = work_dir / "soffice-profile"
    profile.mkdir()
    subprocess.run(
        [
            find_soffice(),
            f"-env:UserInstallation={profile.resolve().as_uri()}",
            "--headless",
            "--convert-to",
            "xlsx",
            "--outdir",
            str(work_dir),
            str(source_xls),
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    converted = work_dir / "GSE49543_further_processed.xlsx"
    if not converted.exists():
        raise FileNotFoundError(converted)
    return converted


def normalize_excel_date(value: object) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, pd.Timestamp):
        return f"{value.day}-{MONTHS[value.month - 1]}"
    text = str(value).strip()
    match = re.fullmatch(r"(?:2013-(\d{2})-(\d{2}))(?: 00:00:00)?", text)
    if match:
        return f"{int(match.group(2))}-{MONTHS[int(match.group(1)) - 1]}"
    return text


def load_probe_rows(source_gz: Path, work_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    workbook = convert_workbook(source_gz, work_dir)
    raw = pd.read_excel(
        workbook,
        sheet_name="Group-YC-MA-MP-SP",
        header=None,
        dtype=object,
    )
    header = raw.iloc[4].tolist()
    symbol_column = next(
        index
        for index, value in enumerate(header)
        if str(value).strip() == "Gene Symbol"
    )
    sample_columns = [
        index
        for index, value in enumerate(raw.iloc[3].tolist())
        if re.fullmatch(r"(?:YC|MA|MP|SP)-\d+", str(value).strip())
    ]
    if len(sample_columns) != 41:
        raise ValueError(
            f"Expected 41 GSE49543 sample columns, found {len(sample_columns)}."
        )
    data = raw.iloc[5:].copy()
    data.columns = range(data.shape[1])
    symbols = data[symbol_column].map(normalize_excel_date).str.strip()
    values = data[sample_columns].apply(pd.to_numeric, errors="coerce")
    valid = (
        values.notna().all(axis=1)
        & symbols.ne("")
        & symbols.ne("nan")
    )
    symbols = symbols[valid].reset_index(drop=True)
    values = values[valid].reset_index(drop=True)
    sample_names = [str(raw.iloc[3, index]).strip() for index in sample_columns]
    values.columns = sample_names
    probe_rows = symbols.to_frame("symbol")
    probe_rows = pd.concat([probe_rows, values], axis=1)
    return probe_rows, pd.DataFrame(
        {
            "source_row": np.flatnonzero(valid.to_numpy()) + 6,
            "symbol": symbols,
            "probe_mean_linear_signal": values.mean(axis=1),
        }
    )


def write_scheme(
    probe_rows: pd.DataFrame,
    output: Path,
    scheme: str,
) -> pd.DataFrame:
    sample_names = [
        column
        for column in probe_rows.columns
        if column not in {"symbol", "source_order"}
    ]
    if scheme == "first_probe":
        selected = probe_rows[~probe_rows["symbol"].duplicated(keep="first")]
        linear = selected.set_index("symbol")[sample_names]
    elif scheme == "mean_probe":
        linear = probe_rows.groupby("symbol", sort=True)[sample_names].mean()
    elif scheme == "max_mean_probe":
        ranked = probe_rows.copy()
        ranked["probe_mean"] = ranked[sample_names].mean(axis=1)
        ranked = ranked.sort_values(
            ["symbol", "probe_mean", "source_order"],
            ascending=[True, False, True],
        )
        selected = ranked[~ranked["symbol"].duplicated(keep="first")]
        linear = selected.set_index("symbol")[sample_names]
    else:
        raise ValueError(scheme)
    matrix = np.log2(linear.astype(float) + 1.0).sort_index()
    matrix.index.name = "Symbol"
    matrix.to_csv(output / f"GSE49543_{scheme}_log2.csv", encoding="utf-8")
    return matrix


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geo-dir", type=Path, default=DEFAULT_GEO_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    geo_dir = args.geo_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = geo_dir / "GSE49543_further_processed.xls.gz"
    if not source.exists():
        raise FileNotFoundError(source)

    with tempfile.TemporaryDirectory(prefix="v14_probe_sensitivity_") as temp:
        probe_rows, probe_audit = load_probe_rows(source, Path(temp))
    probe_rows["source_order"] = np.arange(len(probe_rows))
    audit_path = output / "GSE49543_probe_audit.csv"
    probe_audit.to_csv(audit_path, index=False, encoding="utf-8")

    matrices = {
        scheme: write_scheme(probe_rows, output, scheme)
        for scheme in ("first_probe", "mean_probe", "max_mean_probe")
    }
    duplicate_symbols = int(probe_rows["symbol"].duplicated().sum())
    manifest = {
        "generator": "code/python/prepare_gse49543_probe_sensitivity.py",
        "source_file": source.name,
        "source_sha256": sha256(source),
        "probe_rows": int(len(probe_rows)),
        "duplicate_symbol_rows": duplicate_symbols,
        "sample_columns": int(probe_rows.shape[1] - 2),
        "rules": {
            scheme: {
                "description": {
                    "first_probe": "First valid source row for each exact gene symbol.",
                    "mean_probe": "Arithmetic mean of linear RMA signals across rows sharing a symbol.",
                    "max_mean_probe": (
                        "Row with the highest mean linear RMA signal across all 41 samples "
                        "for each symbol; ties use source order."
                    ),
                }[scheme],
                "rows": int(matrices[scheme].shape[0]),
                "output": f"GSE49543_{scheme}_log2.csv",
                "output_sha256": sha256(output / f"GSE49543_{scheme}_log2.csv"),
            }
            for scheme in matrices
        },
    }
    (output / "GSE49543_probe_sensitivity_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
