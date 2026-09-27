#!/usr/bin/env python3
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import (
    GridSearchCV,
    RepeatedStratifiedKFold,
    StratifiedKFold,
    cross_val_predict,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


ROOT = Path("/Users/jijiabiao/Desktop/耳聋课题")
OUT = ROOT / "2026_supplement_analysis"
TEMP = Path("/tmp/ear_study_learning")
EXPRESSION = ROOT / "2022耳聋耳鸣" / "public" / "exp.txt"
warnings.filterwarnings("ignore")


class DatasetZScore(BaseEstimator, TransformerMixin):
    def fit(self, x, y=None):
        x = np.asarray(x, dtype=float)
        self.mean_ = x.mean(axis=0)
        self.scale_ = x.std(axis=0)
        self.scale_[self.scale_ == 0] = 1.0
        return self

    def transform(self, x):
        x = np.asarray(x, dtype=float)
        return (x - self.mean_) / self.scale_


def build_pipeline(model_name: str, k: int) -> Pipeline:
    if model_name == "logistic":
        model = LogisticRegression(
            C=1.0,
            penalty="l2",
            solver="liblinear",
            class_weight="balanced",
            max_iter=5000,
            random_state=2026,
        )
    elif model_name == "svm":
        model = SVC(
            C=1.0,
            kernel="linear",
            class_weight="balanced",
            probability=True,
            random_state=2026,
        )
    else:
        raise ValueError(model_name)
    return Pipeline(
        [
            ("zscore", DatasetZScore()),
            ("select", SelectKBest(score_func=f_classif, k=k)),
            ("model", model),
        ]
    )


def nested_binary_cv(
    x: np.ndarray,
    y: np.ndarray,
    name: str,
) -> tuple[list[dict], dict[str, object]]:
    outer = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=2026)
    inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=2026)
    grid = {
        "select__k": [10, 20, 50, 100],
        "model__C": [0.1, 1.0, 10.0] if True else [],
    }
    rows = []
    model_selections = []
    fold_index = 0
    for train_index, test_index in outer.split(x, y):
        fold_index += 1
        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        for model_name in ("logistic", "svm"):
            estimator = build_pipeline(model_name, 20)
            search = GridSearchCV(
                estimator,
                grid,
                scoring="roc_auc",
                cv=inner,
                n_jobs=-1,
                refit=True,
            )
            search.fit(x_train, y_train)
            probabilities = search.predict_proba(x_test)[:, 1]
            predictions = search.predict(x_test)
            rows.append(
                {
                    "task": name,
                    "fold": fold_index,
                    "model": model_name,
                    "auc": roc_auc_score(y_test, probabilities),
                    "accuracy": accuracy_score(y_test, predictions),
                    "f1": f1_score(y_test, predictions),
                    "best_params": search.best_params_,
                }
            )
            selected = search.best_estimator_.named_steps["select"].get_support()
            model_selections.append(
                {
                    "task": name,
                    "fold": fold_index,
                    "model": model_name,
                    "selected_indices": np.flatnonzero(selected).tolist(),
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(
        OUT / f"nested_cv_{name}_folds.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary = {
        "task": name,
        "models": {},
    }
    for model_name, group in result.groupby("model"):
        summary["models"][model_name] = {
            "auc_mean": float(group["auc"].mean()),
            "auc_sd": float(group["auc"].std(ddof=1)),
            "accuracy_mean": float(group["accuracy"].mean()),
            "f1_mean": float(group["f1"].mean()),
        }
    return model_selections, summary


def external_validation(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    y_test: np.ndarray,
) -> dict[str, float]:
    results = {}
    for model_name in ("logistic", "svm"):
        estimator = build_pipeline(model_name, 20)
        estimator.fit(x_train, y_train)
        probabilities = estimator.predict_proba(x_test)[:, 1]
        predictions = estimator.predict(x_test)
        results[model_name] = {
            "auc": float(roc_auc_score(y_test, probabilities)),
            "accuracy": float(accuracy_score(y_test, predictions)),
            "f1": float(f1_score(y_test, predictions)),
        }
    return results


def read_gse49543() -> tuple[pd.DataFrame, pd.Series]:
    expression = pd.read_csv(EXPRESSION, sep="\t", index_col=0)
    expression = expression.groupby(level=0).mean()
    logged = np.log2(expression + 1.0)
    group = pd.Series(
        [column.split("-")[0] for column in logged.columns],
        index=logged.columns,
    )
    return logged, group


def read_gse233798() -> tuple[pd.DataFrame, pd.Series]:
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
    expression = np.log2(frame[old + young] + 1.0).T
    labels = pd.Series(["old"] * len(old) + ["young"] * len(young), index=expression.index)
    expression.columns = expression.columns.astype(str)
    return expression, labels


def main() -> int:
    expression, group = read_gse49543()
    tasks = {
        "severe_vs_rest": (group == "SP").astype(int),
        "hearing_loss_vs_normal": group.isin(["MP", "SP"]).astype(int),
    }
    summaries = {}
    selections = {}
    for task_name, labels in tasks.items():
        selections[task_name], summaries[task_name] = nested_binary_cv(
            expression.to_numpy(dtype=float).T,
            labels.to_numpy(dtype=int),
            task_name,
        )

    external_expression, external_group = read_gse233798()
    common_genes = expression.index.intersection(external_expression.columns)
    x_train = expression.loc[common_genes].to_numpy(dtype=float).T
    x_test = external_expression[common_genes].to_numpy(dtype=float)
    y_train = group.isin(["MP", "SP"]).astype(int).to_numpy()
    y_test = (external_group == "old").astype(int).to_numpy()
    external = external_validation(x_train, y_train, x_test, y_test)

    selected_counts = {}
    gene_names = expression.index.to_numpy()
    for task_name, entries in selections.items():
        for model_name in ("logistic", "svm"):
            counts: dict[str, int] = {}
            for entry in entries:
                if entry["model"] != model_name:
                    continue
                for index in entry["selected_indices"]:
                    gene = str(gene_names[index])
                    counts[gene] = counts.get(gene, 0) + 1
            selected_counts[f"{task_name}_{model_name}"] = sorted(
                counts.items(),
                key=lambda item: (-item[1], item[0]),
            )[:100]

    summary = {
        "nested_cv": summaries,
        "external_GSE233798": external,
        "external_common_genes": len(common_genes),
        "selected_feature_stability": selected_counts,
    }
    (OUT / "nested_cv_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
