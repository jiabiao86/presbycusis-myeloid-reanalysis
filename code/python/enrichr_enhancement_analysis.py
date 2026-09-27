#!/usr/bin/env python3
"""Run TF, miRNA, drug, and second-pass functional enrichment with Enrichr."""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd


ROOT = Path("/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis")
OUT = ROOT / "enhancement_analysis"
ENRICHR_ADD = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH = "https://maayanlab.cloud/Enrichr/enrich"

GENE_SET_FILES = {
    "M12_module": ROOT / "wgcna_module_assignment.csv",
    "Severity_positive_FDR005": ROOT / "limma_severity_trend.csv",
}

LIBRARIES = {
    "transcription_factors": "ChEA_2022",
    "micrornas": "TargetScan_microRNA_2017",
    "drug_repositioning": "DSigDB",
    "go_biological_process": "GO_Biological_Process_2025",
    "kegg_mouse": "KEGG_2019_Mouse",
}


def derive_gene_sets() -> dict[str, list[str]]:
    modules = pd.read_csv(GENE_SET_FILES["M12_module"], encoding="utf-8-sig")
    m12 = modules.loc[modules["module"] == "M12", "gene"].dropna().astype(str).tolist()

    trend = pd.read_csv(GENE_SET_FILES["Severity_positive_FDR005"], encoding="utf-8-sig")
    positive = trend.loc[
        (trend["adj.P.Val"] < 0.05) & (trend["logFC"] > 0),
        "gene",
    ].dropna().astype(str).tolist()

    clean = lambda genes: sorted(
        {
            gene.upper()
            for value in genes
            for gene in value.split(" /// ")
            if gene and gene != "---"
        }
    )
    return {"M12_module": clean(m12), "Severity_positive_FDR005": clean(positive)}


def request_json(url: str, payload: bytes | None = None, content_type: str | None = None):
    headers = {"User-Agent": "Codex academic analysis"}
    if content_type:
        headers["Content-Type"] = content_type
    request = urllib.request.Request(
        url,
        data=payload,
        headers=headers,
    )
    with urllib.request.urlopen(request, timeout=90) as response:
        return json.loads(response.read().decode("utf-8"))


def enrich(gene_list: list[str], description: str, library: str) -> pd.DataFrame:
    boundary = "----CodexEnrichrBoundary"
    fields = {
        "list": "\n".join(gene_list),
        "description": description,
    }
    body_parts = []
    for name, value in fields.items():
        body_parts.append(f"--{boundary}\r\n".encode())
        body_parts.append(f'Content-Disposition: form-data; name="{name}"\r\n\r\n'.encode())
        body_parts.append(str(value).encode())
        body_parts.append(b"\r\n")
    body_parts.append(f"--{boundary}--\r\n".encode())
    payload = b"".join(body_parts)
    add_result = request_json(
        ENRICHR_ADD,
        payload,
        content_type=f"multipart/form-data; boundary={boundary}",
    )
    user_list_id = add_result["userListId"]
    time.sleep(0.4)
    query = urllib.parse.urlencode(
        {"userListId": user_list_id, "backgroundType": library}
    )
    result = request_json(f"{ENRICHR_ENRICH}?{query}")
    terms = result.get(library, [])
    frame = pd.DataFrame(
        [
            {
                "rank": rank,
                "term": row[1],
                "p_value": row[2],
                "z_score": row[3],
                "combined_score": row[4],
                "overlapping_genes": ";".join(row[5]),
                "adjusted_p_value": row[6],
                "old_p_value": row[7],
                "old_adjusted_p_value": row[8],
            }
            for rank, row in enumerate(terms, start=1)
        ]
    )
    return frame


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gene_sets = derive_gene_sets()
    summary: dict[str, dict[str, object]] = {}

    for set_name, genes in gene_sets.items():
        summary[set_name] = {"gene_count": len(genes), "libraries": {}}
        for analysis_name, library in LIBRARIES.items():
            frame = enrich(genes, set_name, library)
            output = OUT / f"{set_name}_{analysis_name}.csv"
            frame.to_csv(output, index=False, encoding="utf-8-sig")
            summary[set_name]["libraries"][analysis_name] = {
                "library": library,
                "terms": int(len(frame)),
                "top_terms": frame.head(10).to_dict(orient="records"),
            }
            print(f"{set_name} / {analysis_name}: {len(frame)} terms")

    (OUT / "enhancement_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True),
        encoding="utf-8",
    )
    print(f"Wrote enrichment results to {OUT}")


if __name__ == "__main__":
    main()
