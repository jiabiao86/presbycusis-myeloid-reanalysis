#!/usr/bin/env python3
"""Build the 80-sample V14 mouse metadata audit from deposited GEO files."""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
GEO = REPO_ROOT / "inputs" / "geo"
PREPARED = Path(
    os.environ.get("V14_PREPARED_DIR", REPO_ROOT / "inputs" / "prepared")
).resolve()
RESULTS = Path(
    os.environ.get("V14_PRIMARY_RESULTS_DIR", REPO_ROOT / "results" / "primary")
).resolve()
RETRIEVAL_DATE = "2026-10-29"

GROUP_NAME = {
    "YC": "Young Control",
    "MA": "Middle-Aged",
    "MP": "Mild Presbycusis",
    "SP": "Severe Presbycusis",
}
GROUP_CODE = {value: key for key, value in GROUP_NAME.items()}
GROUP_CODE["Sever Presbycusis"] = "SP"

AGE_SUMMARY = {
    "YC": (3.5, 0.4),
    "MA": (12.3, 1.3),
    "MP": (27.7, 3.4),
    "SP": (30.6, 1.9),
}

SERIES = {
    "GSE49543": {
        "tissue": "cochlea",
        "matrix": GEO / "GSE49543_series_matrix.txt.gz",
        "expression": PREPARED / "expression_GSE49543_log2.csv",
        "documented_n": 40,
    },
    "GSE49522": {
        "tissue": "inferior_colliculus",
        "matrix": GEO / "GSE49522_series_matrix.txt.gz",
        "expression": PREPARED / "expression_GSE49522_log2.csv",
        "documented_n": 39,
    },
}


def strip_quotes(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] == '"':
        return value[1:-1]
    return value


def read_series_matrix(path: Path) -> tuple[list[dict[str, str]], dict[str, str]]:
    wanted = {
        "!Sample_title": "title",
        "!Sample_geo_accession": "gsm",
    }
    character_rows: list[list[str]] = []
    parsed: dict[str, list[str]] = {}
    series_meta: dict[str, str] = {}

    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!"):
                continue
            fields = [strip_quotes(value) for value in line.split("\t")]
            key = fields[0]
            if key in wanted:
                parsed[wanted[key]] = fields[1:]
            elif key == "!Sample_characteristics_ch1":
                character_rows.append(fields[1:])
            elif key in {"!Series_title", "!Series_overall_design"}:
                series_meta[key] = fields[1] if len(fields) > 1 else "NA"

    if not parsed.get("title") or not parsed.get("gsm"):
        raise ValueError(f"Could not parse sample titles or accessions from {path}")
    if len(parsed["title"]) != len(parsed["gsm"]):
        raise ValueError(f"Sample title/accession length mismatch in {path}")
    for row in character_rows:
        if len(row) != len(parsed["title"]):
            raise ValueError(f"Sample characteristics length mismatch in {path}")

    records: list[dict[str, str]] = []
    for index, (title, gsm) in enumerate(zip(parsed["title"], parsed["gsm"])):
        characteristics: dict[str, str] = {}
        for row in character_rows:
            item = row[index]
            if ": " in item:
                key, value = item.split(": ", 1)
                characteristics[key.lower()] = value
        match = re.search(r"Mice-(\d+)$", title)
        if not match:
            raise ValueError(f"Could not extract animal ID from sample title: {title}")
        records.append(
            {
                "title": title,
                "gsm": gsm,
                "animal_id": match.group(1),
                "reported_group": characteristics.get("hearing status", "NA"),
                "sex": characteristics.get("gender", "NA"),
            }
        )
    return records, series_meta


def expression_columns(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
    return header[1:]


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    all_records: list[dict[str, str]] = []
    audit: dict[str, object] = {
        "retrieval_date": RETRIEVAL_DATE,
        "datasets": {},
        "deposited_animal_level_fields": {
            "age_months": False,
            "abr_threshold_db": False,
            "dpoae_db": False,
        },
        "notes": [
            "GEO further-processed workbooks contain expression values and cohort labels only.",
            "No per-animal numeric age, ABR threshold, or DPOAE value was found.",
            "The GSE49543 overall design states that each mouse had DPOAE amplitude, but no per-animal DPOAE values are present in the deposited series matrix, family SOFT record, or further-processed workbook.",
            "Published cohort age means and SDs are retained only as cohort metadata.",
        ],
    }

    for accession, config in SERIES.items():
        records, series_meta = read_series_matrix(config["matrix"])
        expr_cols = expression_columns(config["expression"])
        expected_expr_cols = {
            f"{GROUP_CODE[record['reported_group']]}-{record['animal_id']}"
            if accession == "GSE49543"
            else record["gsm"]
            for record in records
        }
        if set(expr_cols) != expected_expr_cols:
            missing = sorted(expected_expr_cols - set(expr_cols))
            extra = sorted(set(expr_cols) - expected_expr_cols)
            raise ValueError(
                f"{accession} expression/GSM mismatch; missing={missing}, extra={extra}"
            )

        group_counts: dict[str, int] = {}
        sex_counts: dict[str, dict[str, int]] = {}
        for record in records:
            report = record["reported_group"]
            if report not in GROUP_CODE:
                raise ValueError(f"Unrecognized group label in {accession}: {report}")
            code = GROUP_CODE[report]
            age_mean, age_sd = AGE_SUMMARY[code]
            record.update(
                {
                    "dataset": accession,
                    "tissue": config["tissue"],
                    "adjudicated_group": GROUP_NAME[code],
                    "cohort_group": code,
                    "age_mean_months": age_mean,
                    "age_sd_months": age_sd,
                    "expression_column": (
                        f"{code}-{record['animal_id']}"
                        if accession == "GSE49543"
                        else record["gsm"]
                    ),
                }
            )
            group_counts[code] = group_counts.get(code, 0) + 1
            sex_counts.setdefault(code, {})
            sex_counts[code][record["sex"]] = sex_counts[code].get(record["sex"], 0) + 1
            all_records.append(record)
        audit["datasets"][accession] = {
            "tissue": config["tissue"],
            "observed_n": len(records),
            "documented_n": config["documented_n"],
            "counts_by_group": dict(sorted(group_counts.items())),
            "sex_by_group": {key: dict(sorted(value.items())) for key, value in sorted(sex_counts.items())},
            "sex_complete": all(record["sex"] != "NA" for record in records),
            "series_title": series_meta.get("!Series_title", "NA"),
            "series_overall_design": series_meta.get("!Series_overall_design", "NA"),
            "count_discrepancy": len(records) != config["documented_n"],
            "count_discrepancy_note": (
                "Observed 41 GSMs versus documented 40 (YC 8, MA 17, MP 9, SP 6); "
                "the two unmatched cochlear records are animal 13 (SP, GSM1201602) "
                "and animal 16 (YC, GSM1201568). Both remain in the primary analysis."
                if accession == "GSE49543"
                else "Observed count matches the documented 39."
            ),
        }

    pair_map: dict[str, dict[str, str]] = {}
    for record in all_records:
        pair_map.setdefault(record["dataset"], {})
        pair_map[record["dataset"]][record["animal_id"]] = record["gsm"]
    cochlea_ids = set(pair_map.get("GSE49543", {}))
    ic_ids = set(pair_map.get("GSE49522", {}))
    shared_ids = cochlea_ids & ic_ids
    audit["cross_tissue"] = {
        "cochlea_animals": len(cochlea_ids),
        "inferior_colliculus_animals": len(ic_ids),
        "paired_animals": len(shared_ids),
        "cochlea_only_animal_ids": sorted(cochlea_ids - ic_ids, key=int),
        "inferior_colliculus_only_animal_ids": sorted(ic_ids - cochlea_ids, key=int),
    }

    fields = [
        "sample_id",
        "dataset",
        "tissue",
        "gsm",
        "animal_id",
        "sample_title",
        "expression_column",
        "reported_group",
        "adjudicated_group",
        "cohort_group",
        "sex",
        "age_months",
        "age_mean_months",
        "age_sd_months",
        "age_scope",
        "abr_threshold_db",
        "dpoae_db",
        "hearing_data_status",
        "source_record_url",
        "source_publication_original",
        "source_publication_reuse_or_related",
        "source_publication_classification",
        "source_retrieved_on",
        "missing_reason",
        "include_main",
        "inclusion_reason",
        "exclusion_reason",
        "design_count_mismatch",
        "count_audit_note",
        "paired_tissue_gsm",
        "paired_tissue_available",
        "notes",
    ]
    rows: list[dict[str, object]] = []
    for record in all_records:
        accession = record["dataset"]
        paired_gsm = pair_map.get(
            "GSE49522" if accession == "GSE49543" else "GSE49543", {}
        ).get(record["animal_id"], "NA")
        count_note = (
            "GSE49543 contains 41 observed GSMs versus 40 documented. The unmatched "
            "cochlear records are animal 13 (SP, GSM1201602) and animal 16 "
            "(YC, GSM1201568); both are retained in the primary analysis."
            if accession == "GSE49543"
            else "GSE49522 contains 39 observed GSMs, matching its documented total."
        )
        notes = []
        if record["reported_group"] == "Sever Presbycusis":
            notes.append("Spelling normalized to Severe Presbycusis in adjudicated_group only.")
        if paired_gsm == "NA":
            notes.append("No same-animal sample was found in the other tissue series.")
            if accession == "GSE49543" and record["animal_id"] == "13":
                notes.append("Unmatched cochlear record is SP according to the source label.")
            if accession == "GSE49543" and record["animal_id"] == "16":
                notes.append("Unmatched cochlear record is YC according to the source label.")
        rows.append(
            {
                "sample_id": record["gsm"],
                "dataset": accession,
                "tissue": record["tissue"],
                "gsm": record["gsm"],
                "animal_id": record["animal_id"],
                "sample_title": record["title"],
                "expression_column": record["expression_column"],
                "reported_group": record["reported_group"],
                "adjudicated_group": record["adjudicated_group"],
                "cohort_group": record["cohort_group"],
                "sex": record["sex"],
                "age_months": "NA",
                "age_mean_months": record["age_mean_months"],
                "age_sd_months": record["age_sd_months"],
                "age_scope": "published_cohort_mean_sd_only",
                "abr_threshold_db": "NA",
                "dpoae_db": "NA",
                "hearing_data_status": "animal_level_not_deposited",
                "source_record_url": (
                    f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={record['gsm']}"
                ),
                "source_publication_original": (
                    "D'Souza et al. (2008); PMID 18455804; "
                    "DOI 10.1016/j.jneumeth.2008.02.022"
                    if accession == "GSE49543"
                    else (
                        "D'Souza et al. (2013) dataset citation; "
                        "NCBI Gene Expression Omnibus GSE49522; "
                        "BioProject PRJNA214236"
                    )
                ),
                "source_publication_reuse_or_related": (
                    "Tadros et al. (2014); PMID 24587312; PMC3938674; "
                    "DOI 10.1371/journal.pone.0090279; "
                    "Pineros et al. (2024); PMID 38359484; PMC10939722; "
                    "DOI 10.1016/j.heares.2024.108971"
                    if accession == "GSE49543"
                    else "None listed in the GEO/BioProject record"
                ),
                "source_publication_classification": (
                    "Original dataset publication identified; later GEO-linked "
                    "publications recorded separately as reuse or related sources"
                    if accession == "GSE49543"
                    else (
                        "Dataset citation identified; no linked journal article "
                        "in GEO/BioProject records"
                    )
                ),
                "source_retrieved_on": RETRIEVAL_DATE,
                "missing_reason": (
                    "age_months: only published cohort mean and SD available; "
                    "abr_threshold_db: per-animal value not deposited; "
                    "dpoae_db: per-animal value not deposited"
                ),
                "include_main": "TRUE",
                "inclusion_reason": (
                    "Observed GSM with complete expression column and deposited group/sex metadata; retained in the primary observed-sample analysis."
                ),
                "exclusion_reason": "NA",
                "design_count_mismatch": "TRUE" if accession == "GSE49543" else "FALSE",
                "count_audit_note": count_note,
                "paired_tissue_gsm": paired_gsm,
                "paired_tissue_available": "TRUE" if paired_gsm != "NA" else "FALSE",
                "notes": " ".join(notes) if notes else "NA",
            }
        )

    PREPARED.mkdir(parents=True, exist_ok=True)
    out_path = PREPARED / "sample_master.csv"
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    cochlea_groups = {
        row["animal_id"]: row["cohort_group"]
        for row in rows
        if row["dataset"] == "GSE49543"
    }
    ic_groups = {
        row["animal_id"]: row["cohort_group"]
        for row in rows
        if row["dataset"] == "GSE49522"
    }
    shared_animal_ids = sorted(
        set(cochlea_groups) & set(ic_groups), key=int
    )
    audit["cross_tissue"]["paired_group_concordant"] = sum(
        cochlea_groups[animal_id] == ic_groups[animal_id]
        for animal_id in shared_animal_ids
    )
    audit["cross_tissue"]["paired_group_discordant"] = sum(
        cochlea_groups[animal_id] != ic_groups[animal_id]
        for animal_id in shared_animal_ids
    )
    (RESULTS / "sample_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=True) + "\n", encoding="utf-8"
    )
    age_fields = [
        "cohort_group",
        "age_mean_months",
        "age_sd_months",
        "mean_minus_one_sd_months",
        "mean_plus_one_sd_months",
        "individual_age_available",
        "comparability_note",
    ]
    age_rows = []
    for code, (mean, sd) in AGE_SUMMARY.items():
        age_rows.append(
            {
                "cohort_group": code,
                "age_mean_months": mean,
                "age_sd_months": sd,
                "mean_minus_one_sd_months": round(mean - sd, 1),
                "mean_plus_one_sd_months": round(mean + sd, 1),
                "individual_age_available": "FALSE",
                "comparability_note": (
                    "MP and SP cohort-level mean +/- 1 SD ranges overlap, but no individual ages were deposited."
                    if code in {"MP", "SP"}
                    else "Published cohort summary only; no individual ages were deposited."
                ),
            }
        )
    with (RESULTS / "age_comparability.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=age_fields)
        writer.writeheader()
        writer.writerows(age_rows)

    if len(rows) != 80:
        raise ValueError(f"Expected 80 sample rows, wrote {len(rows)}")
    print(f"Wrote {out_path.relative_to(REPO_ROOT)} with {len(rows)} samples")
    print(json.dumps(audit, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
