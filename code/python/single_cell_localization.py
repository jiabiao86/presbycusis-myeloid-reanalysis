#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu
from sklearn.cluster import KMeans


ROOT = Path("/Users/jijiabiao/Desktop/耳聋课题")
OUT = ROOT / "2026_supplement_analysis"
SC_DIR = Path("/tmp/ear_study_learning/GSE274279")
WORK = Path("/tmp/ear_study_learning/GSE274279_10x")


CANDIDATES = [
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
    "C1qc",
    "C3ar1",
    "Ms4a7",
    "Apoe",
    "Mif",
    "Ddx1",
    "Lamp3",
    "Tbk1",
    "Fgf16",
    "Csf1r",
    "P2ry12",
    "Trem2",
]

MARKERS = {
    "Hair cells": ["Pou4f3", "Myo7a", "Atoh1", "Pcp4", "Gfi1"],
    "Inner hair cells": ["Slc17a8", "Ocm", "Calb2"],
    "Outer hair cells": ["Slc26a5", "Kcnq4", "Cdh23"],
    "Supporting cells": ["Gjb2", "Sox2", "Fbxo2", "Slc1a3", "S100a6"],
    "Spiral ganglion neurons": ["Nefh", "Snap25", "Tubb3", "Nefl"],
    "Glia/Schwann": ["Mog", "Mpz", "Plp1", "Sox10"],
    "Stria vascularis": ["Kcnq1", "Slc12a2", "Atp1a1", "Dct", "Ednrb"],
    "Fibrocytes": ["Slc26a4", "Col1a1", "Col1a2", "Dcn"],
    "Endothelial": ["Pecam1", "Cldn5", "Vwf"],
    "Macrophages/Microglia": [
        "Ptprc",
        "Cd68",
        "Cd74",
        "C1qa",
        "C1qb",
        "Tyrobp",
        "Laptm5",
        "Fcgr3",
        "Ctss",
        "Ms4a7",
        "Mpeg1",
        "Csf1r",
        "P2ry12",
        "Trem2",
    ],
    "T cells": ["Cd3d", "Cd3e", "Cd3g", "Lck"],
    "B cells": ["Cd79a", "Cd79b", "Ms4a1"],
}


def prepare_10x() -> list[tuple[str, int, Path]]:
    WORK.mkdir(parents=True, exist_ok=True)
    samples = [
        ("3M", 3, SC_DIR / "3M_matrix.mtx.gz", SC_DIR / "3M_barcodes.tsv.gz"),
        ("12M", 12, SC_DIR / "12M_matrix.mtx.gz", SC_DIR / "12M_barcodes.tsv.gz"),
        ("24M", 24, SC_DIR / "24M_matrix.mtx.gz", SC_DIR / "24M_barcodes.tsv.gz"),
    ]
    result = []
    for name, age, matrix, barcodes in samples:
        sample_dir = WORK / name
        sample_dir.mkdir(exist_ok=True)
        links = {
            "matrix.mtx.gz": matrix,
            "barcodes.tsv.gz": barcodes,
            "features.tsv.gz": SC_DIR / "features.tsv.gz",
        }
        for target_name, source in links.items():
            target = sample_dir / target_name
            if target.exists() or target.is_symlink():
                target.unlink()
            os.symlink(source, target)
        result.append((name, age, sample_dir))
    return result


def bh(p_values: np.ndarray) -> np.ndarray:
    p_values = np.clip(np.asarray(p_values, dtype=float), np.nextafter(0.0, 1.0), 1.0)
    count = len(p_values)
    order = np.argsort(p_values)
    adjusted = p_values[order] * count / (np.arange(count) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(count, dtype=float)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def expression_values(adata, gene: str, mask: np.ndarray) -> np.ndarray:
    if gene not in adata.var_names:
        return np.array([])
    values = adata[:, gene].X
    if hasattr(values, "toarray"):
        values = values.toarray()
    return np.asarray(values, dtype=float).reshape(-1)[mask]


def main() -> int:
    samples = prepare_10x()
    datasets = []
    for name, age, sample_dir in samples:
        adata = sc.read_10x_mtx(
            sample_dir,
            var_names="gene_symbols",
            make_unique=True,
        )
        adata.obs["sample"] = name
        adata.obs["age_months"] = age
        datasets.append(adata)
    adata = sc.concat(datasets, join="inner", label="batch", index_unique="-")
    adata.var_names_make_unique()

    adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[
        (adata.obs["pct_counts_mt"] < 20)
        & (adata.obs["n_genes_by_counts"] < 8000)
    ].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(adata, n_top_genes=2500, batch_key="batch")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.combat(adata_hvg, key="batch")
    adata_hvg.X = adata_hvg.X.copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=40, random_state=2026)

    kmeans = KMeans(n_clusters=16, n_init=50, random_state=2026)
    adata.obs["cluster"] = pd.Categorical(
        kmeans.fit_predict(adata_hvg.obsm["X_pca"][:, :30]).astype(str)
    )
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")
    sc.tl.umap(adata, random_state=2026)

    score_rows = []
    for cell_type, markers in MARKERS.items():
        available = [gene for gene in markers if gene in adata.var_names]
        if not available:
            continue
        sc.tl.score_genes(adata, available, score_name=f"score_{cell_type}")
        for cluster, group in adata.obs.groupby("cluster", observed=True):
            score_rows.append(
                {
                    "cluster": str(cluster),
                    "cell_type": cell_type,
                    "score": float(group[f"score_{cell_type}"].mean()),
                    "markers": len(available),
                }
            )
    score_frame = pd.DataFrame(score_rows)
    score_frame.to_csv(
        OUT / "GSE274279_cluster_marker_scores.csv",
        index=False,
        encoding="utf-8-sig",
    )
    best = (
        score_frame.sort_values(["cluster", "score"], ascending=[True, False])
        .drop_duplicates("cluster")
        .set_index("cluster")["cell_type"]
        .to_dict()
    )
    adata.obs["cell_type"] = adata.obs["cluster"].astype(str).map(best)
    cluster_summary = (
        adata.obs.groupby(["cluster", "cell_type", "age_months"], observed=True)
        .size()
        .rename("cells")
        .reset_index()
    )
    cluster_summary.to_csv(
        OUT / "GSE274279_cell_type_counts.csv",
        index=False,
        encoding="utf-8-sig",
    )

    expression_rows = []
    for gene in CANDIDATES:
        if gene not in adata.var_names:
            continue
        for (cell_type, age), group in adata.obs.groupby(
            ["cell_type", "age_months"],
            observed=True,
        ):
            mask = group.index.to_numpy()
            values = expression_values(adata, gene, adata.obs.index.isin(mask))
            if len(values) == 0:
                continue
            expression_rows.append(
                {
                    "gene": gene,
                    "cell_type": cell_type,
                    "age_months": int(age),
                    "cells": len(values),
                    "mean_log_expression": float(values.mean()),
                    "percent_expressing": float((values > 0).mean() * 100),
                }
            )
    expression_frame = pd.DataFrame(expression_rows)
    expression_frame.to_csv(
        OUT / "GSE274279_candidate_expression_by_cell_type.csv",
        index=False,
        encoding="utf-8-sig",
    )

    immune_mask = adata.obs["cell_type"].eq("Macrophages/Microglia").to_numpy()
    immune = adata[immune_mask].copy()
    immune_rows = []
    for gene in CANDIDATES:
        if gene not in immune.var_names:
            continue
        young = expression_values(
            immune,
            gene,
            immune.obs["age_months"].eq(3).to_numpy(),
        )
        old = expression_values(
            immune,
            gene,
            immune.obs["age_months"].eq(24).to_numpy(),
        )
        if len(young) < 5 or len(old) < 5:
            continue
        statistic, p_value = mannwhitneyu(old, young, alternative="two-sided")
        immune_rows.append(
            {
                "gene": gene,
                "young_cells": len(young),
                "old_cells": len(old),
                "young_mean": float(young.mean()),
                "old_mean": float(old.mean()),
                "log2FC": float(np.log2((old.mean() + 1e-6) / (young.mean() + 1e-6))),
                "mannwhitney_u": float(statistic),
                "p_value": float(p_value),
            }
        )
    immune_result = pd.DataFrame(immune_rows)
    if not immune_result.empty:
        immune_result["fdr_bh"] = bh(immune_result["p_value"].to_numpy())
        immune_result.sort_values("fdr_bh").to_csv(
            OUT / "GSE274279_macrophage_age_effects.csv",
            index=False,
            encoding="utf-8-sig",
        )

    coordinates = pd.DataFrame(
        adata.obsm["X_umap"],
        columns=["UMAP1", "UMAP2"],
        index=adata.obs_names,
    )
    coordinates["cell_type"] = adata.obs["cell_type"].to_numpy()
    coordinates["age_months"] = adata.obs["age_months"].to_numpy()
    coordinates.to_csv(
        OUT / "GSE274279_umap_coordinates.csv",
        encoding="utf-8-sig",
    )
    summary = {
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "samples": adata.obs.groupby("sample", observed=True).size().to_dict(),
        "cell_types": adata.obs["cell_type"].value_counts().to_dict(),
    }
    (OUT / "GSE274279_single_cell_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
