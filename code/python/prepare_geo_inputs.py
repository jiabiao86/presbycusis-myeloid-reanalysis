#!/usr/bin/env python3
"""Rebuild the frozen GSE49543 and GSE49522 gene-level matrices from GEO files."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
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
GEO_DIR = Path(os.environ.get("V14_GEO_DIR", REPO_ROOT / "inputs" / "geo")).resolve()
DEFAULT_OUTPUT = Path(
    os.environ.get("V14_PREPARED_DIR", REPO_ROOT / "inputs" / "prepared")
).resolve()
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
        "LibreOffice is required to convert the GSE49543 .xls workbook. "
        "Install LibreOffice or set V14_SOFFICE to the soffice executable."
    )


def convert_gse49543(source_gz: Path, work_dir: Path) -> Path:
    source_xls = work_dir / "GSE49543_further_processed.xls"
    with gzip.open(source_gz, "rb") as source, source_xls.open("wb") as target:
        shutil.copyfileobj(source, target)

    profile = work_dir / "soffice-profile"
    profile.mkdir()
    command = [
        find_soffice(),
        f"-env:UserInstallation={profile.resolve().as_uri()}",
        "--headless",
        "--convert-to",
        "xlsx",
        "--outdir",
        str(work_dir),
        str(source_xls),
    ]
    subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    converted = work_dir / "GSE49543_further_processed.xlsx"
    if not converted.exists():
        raise FileNotFoundError(f"LibreOffice did not create {converted}")
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


def prepare_gse49543(source_gz: Path, work_dir: Path) -> pd.DataFrame:
    workbook = convert_gse49543(source_gz, work_dir)
    raw = pd.read_excel(
        workbook,
        sheet_name="Group-YC-MA-MP-SP",
        header=None,
        dtype=object,
    )
    if raw.shape[0] < 6 or raw.shape[1] < 60:
        raise ValueError("Unexpected GSE49543 workbook shape.")

    header = raw.iloc[4].tolist()
    if not any(str(value).strip() == "Gene Symbol" for value in header):
        raise ValueError("GSE49543 Gene Symbol column was not found.")
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
    symbols = data[symbol_column].map(normalize_excel_date)
    values = data[sample_columns].apply(pd.to_numeric, errors="coerce")
    keep = values.notna().all(axis=1) & symbols.ne("") & symbols.ne("nan")
    first_rows = ~symbols[keep].duplicated(keep="first")
    symbols = symbols[keep][first_rows].to_numpy()
    values = values[keep][first_rows].to_numpy(dtype=float)
    sample_names = [str(raw.iloc[3, index]).strip() for index in sample_columns]
    matrix = pd.DataFrame(
        np.log2(values + 1.0),
        index=pd.Index(symbols, name="Symbol"),
        columns=sample_names,
    ).sort_index()
    if matrix.shape != (13649, 41):
        raise ValueError(f"Unexpected GSE49543 prepared shape: {matrix.shape}")
    return matrix


def read_gzip_text(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def prepare_gse49522(family_soft: Path, series_matrix: Path) -> pd.DataFrame:
    family_text = read_gzip_text(family_soft)
    start_marker = "!platform_table_begin"
    end_marker = "!platform_table_end"
    start = family_text.index(start_marker) + len(start_marker)
    end = family_text.index(end_marker, start)
    platform = pd.read_csv(
        io.StringIO(family_text[start:end].strip()),
        sep="\t",
        dtype=str,
    )
    platform["symbol"] = platform["Gene Symbol"].fillna("").str.strip()
    platform = platform[
        platform["symbol"].ne("") & ~platform["symbol"].str.contains("///")
    ][["ID", "symbol"]]

    with gzip.open(series_matrix, "rt", encoding="utf-8-sig") as handle:
        lines = handle.readlines()
    matrix_start = next(
        index
        for index, line in enumerate(lines)
        if line.startswith("!series_matrix_table_begin")
    )
    series = pd.read_csv(
        io.StringIO("".join(lines[matrix_start + 1 :])),
        sep="\t",
        quotechar='"',
    ).dropna(how="all")
    series = series[series["ID_REF"].notna()].set_index("ID_REF")
    series = series.apply(pd.to_numeric, errors="coerce")

    joined = platform.merge(
        series.reset_index(),
        left_on="ID",
        right_on="ID_REF",
        how="inner",
    )
    matrix = np.log2(
        joined.groupby("symbol", sort=True)[series.columns].mean() + 1.0
    )
    matrix.index.name = "Symbol"
    if matrix.shape != (13049, 39):
        raise ValueError(f"Unexpected GSE49522 prepared shape: {matrix.shape}")
    return matrix


def combined_action(source: Path) -> dict[str, str]:
    return {
        "linear_signal_probe_selection": (
            "GSE49543: first probe per exact source Gene Symbol; "
            "GSE49522: drop probes containing '///', then average within symbol"
        ),
        "transform": "log2(linear RMA signal or within-symbol mean signal + 1)",
        "source_file": source.name,
        "source_sha256": sha256(source),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Directory for the two prepared expression CSVs and manifest.",
    )
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    gse49543_source = GEO_DIR / "GSE49543_further_processed.xls.gz"
    gse49522_family = GEO_DIR / "GSE49522_family_soft.txt"
    gse49522_matrix = GEO_DIR / "GSE49522_series_matrix.txt.gz"
    for required in (gse49543_source, gse49522_family, gse49522_matrix):
        if not required.exists():
            raise FileNotFoundError(f"Required GEO source not found: {required}")

    with tempfile.TemporaryDirectory(prefix="v14_geo_prepare_") as temp:
        cochlea = prepare_gse49543(gse49543_source, Path(temp))
    inferior_colliculus = prepare_gse49522(gse49522_family, gse49522_matrix)

    cochlea_path = output / "expression_GSE49543_log2.csv"
    inferior_colliculus_path = output / "expression_GSE49522_log2.csv"
    cochlea.to_csv(cochlea_path, encoding="utf-8")
    inferior_colliculus.to_csv(inferior_colliculus_path, encoding="utf-8")

    manifest = {
        "generator": "code/python/prepare_geo_inputs.py",
        "gse49543": {
            **combined_action(gse49543_source),
            "output": cochlea_path.name,
            "rows": int(cochlea.shape[0]),
            "samples": int(cochlea.shape[1]),
            "output_sha256": sha256(cochlea_path),
        },
        "gse49522": {
            **combined_action(gse49522_family),
            "series_matrix": gse49522_matrix.name,
            "series_matrix_sha256": sha256(gse49522_matrix),
            "output": inferior_colliculus_path.name,
            "rows": int(inferior_colliculus.shape[0]),
            "samples": int(inferior_colliculus.shape[1]),
            "output_sha256": sha256(inferior_colliculus_path),
        },
    }
    (output / "geo_preparation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
