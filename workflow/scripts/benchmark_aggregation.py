import os
import re
import pandas as pd
import numpy as np


def to_list(x):
    if isinstance(x, str):
        return [x]
    return list(x)


native_files = to_list(snakemake.input.native)
fair_files = to_list(snakemake.input.fair)
called_sites_files = to_list(snakemake.input.called_sites)
truth_set_path = snakemake.input.truth_set

out_summary = snakemake.output.summary
out_overall = snakemake.output.overall
out_by_comparison = snakemake.output.by_comparison
out_by_negative_type = snakemake.output.by_negative_type
out_by_tool = snakemake.output.by_tool
out_best_scores = snakemake.output.best_scores
out_called_sites_comp = snakemake.output.called_sites_comp
out_called_sites_sum = snakemake.output.called_sites_sum


def tool_comparison_from_path(path, mode):
    parts = path.replace("\\", "/").split("/")
    for i, part in enumerate(parts):
        if part == mode and i + 2 < len(parts):
            tool = parts[i + 1]
            comparison = parts[i + 2]
            return tool, comparison
    return None, None


def comparison_from_called_sites_path(path):
    parts = path.replace("\\", "/").split("/")
    for i, part in enumerate(parts):
        if part == "coverage" and i + 1 < len(parts):
            return parts[i + 1]
    return None


BEST_METRICS_COLS = [
    "score_column", "is_pvalue", "transform", "auroc", "prauc",
    "best_threshold", "best_threshold_original", "f1", "precision", "recall",
]

ACCURACY_SUMMARY_COLS = [
    "tool", "modification_type", "window", "precision", "recall", "f1",
    "tp", "fp", "fn", "tn", "specificity", "mcc", "auprc", "auroc",
    "called_sites", "total_truth", "total_predicted", "total_negative",
]

ACCURACY_OVERALL_COLS = [
    "tool", "window", "precision", "recall", "f1",
    "tp", "fp", "fn", "tn", "specificity", "mcc", "auprc", "auroc",
    "called_sites", "total_truth", "total_predicted", "total_negative",
]

ACCURACY_BY_COMPARISON_COLS = [
    "tool", "modification_type", "comparison", "window", "precision", "recall", "f1",
    "tp", "fp", "fn", "tn", "specificity", "mcc", "auprc", "auroc",
    "called_sites", "total_truth", "total_negative",
]

BY_NEGATIVE_TYPE_COLS = [
    "tool", "modification_type", "negative_type", "window", "precision", "recall", "f1",
    "tp", "fp", "fn", "tn", "specificity", "auprc", "auroc",
    "called_sites", "total_truth", "total_negative", "total_predicted",
]

BY_TOOL_COLS = [
    "tool", "best_score", "mean_auroc", "std_auroc", "mean_prauc", "std_prauc",
    "mean_f1", "std_f1", "n_comparisons",
]

BEST_SCORES_COLS = [
    "tool", "best_score", "n_comparisons_best",
]


def load_best_metrics(files, mode):
    records = []
    for path in files:
        tool, comparison = tool_comparison_from_path(path, mode)
        if tool is None or comparison is None:
            continue
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            continue
        if df.empty:
            continue
        row = df.iloc[0]
        record = {
            "tool": tool,
            "comparison": comparison,
            "score_column": row.get("score_column", np.nan),
            "is_pvalue": row.get("is_pvalue", False),
            "transform": row.get("transform", "none"),
            "auroc": pd.to_numeric(row.get("auroc", np.nan), errors="coerce"),
            "prauc": pd.to_numeric(row.get("prauc", np.nan), errors="coerce"),
            "best_threshold": pd.to_numeric(row.get("best_threshold", np.nan), errors="coerce"),
            "best_threshold_original": pd.to_numeric(row.get("best_threshold_original", np.nan), errors="coerce"),
            "f1": pd.to_numeric(row.get("f1", np.nan), errors="coerce"),
            "precision": pd.to_numeric(row.get("precision", np.nan), errors="coerce"),
            "recall": pd.to_numeric(row.get("recall", np.nan), errors="coerce"),
        }
        records.append(record)
    return pd.DataFrame(records)


truth = pd.read_csv(truth_set_path, sep="\t")
if "label" in truth.columns:
    truth_pos = truth[truth["label"] != "-"]
else:
    truth_pos = truth
total_truth = len(truth_pos)

native_df = load_best_metrics(native_files, "native")
fair_df = load_best_metrics(fair_files, "fair")

if not native_df.empty:
    primary_df = native_df.copy()
elif not fair_df.empty:
    primary_df = fair_df.copy()
else:
    primary_df = pd.DataFrame()

os.makedirs(os.path.dirname(out_summary), exist_ok=True)

if primary_df.empty:
    pd.DataFrame(columns=ACCURACY_SUMMARY_COLS).to_csv(out_summary, sep="\t", index=False)
    pd.DataFrame(columns=ACCURACY_OVERALL_COLS).to_csv(out_overall, sep="\t", index=False)
    pd.DataFrame(columns=ACCURACY_BY_COMPARISON_COLS).to_csv(out_by_comparison, sep="\t", index=False)
    pd.DataFrame(columns=BY_NEGATIVE_TYPE_COLS).to_csv(out_by_negative_type, sep="\t", index=False)
    pd.DataFrame(columns=BY_TOOL_COLS).to_csv(out_by_tool, sep="\t", index=False)
    pd.DataFrame(columns=BEST_SCORES_COLS).to_csv(out_best_scores, sep="\t", index=False)
    pd.DataFrame(columns=["tool", "comparison", "n_covered"]).to_csv(out_called_sites_comp, sep="\t", index=False)
    pd.DataFrame(columns=["tool", "total_covered", "n_comparisons"]).to_csv(out_called_sites_sum, sep="\t", index=False)
else:
    summary_records = []
    for _, row in primary_df.iterrows():
        rec = {
            "tool": row["tool"],
            "modification_type": "all",
            "window": 0,
            "precision": row["precision"],
            "recall": row["recall"],
            "f1": row["f1"],
            "tp": 0,
            "fp": 0,
            "fn": 0,
            "tn": 0,
            "specificity": np.nan,
            "mcc": np.nan,
            "auprc": row["prauc"],
            "auroc": row["auroc"],
            "called_sites": np.nan,
            "total_truth": total_truth,
            "total_predicted": np.nan,
            "total_negative": np.nan,
        }
        summary_records.append(rec)

    summary_df = pd.DataFrame(summary_records, columns=ACCURACY_SUMMARY_COLS)
    summary_df.to_csv(out_summary, sep="\t", index=False)

    overall_records = []
    for tool, grp in primary_df.groupby("tool"):
        rec = {
            "tool": tool,
            "window": 0,
            "precision": grp["precision"].mean(),
            "recall": grp["recall"].mean(),
            "f1": grp["f1"].mean(),
            "tp": 0,
            "fp": 0,
            "fn": 0,
            "tn": 0,
            "specificity": np.nan,
            "mcc": np.nan,
            "auprc": grp["prauc"].mean(),
            "auroc": grp["auroc"].mean(),
            "called_sites": np.nan,
            "total_truth": total_truth,
            "total_predicted": np.nan,
            "total_negative": np.nan,
        }
        overall_records.append(rec)

    overall_df = pd.DataFrame(overall_records, columns=ACCURACY_OVERALL_COLS)
    overall_df.to_csv(out_overall, sep="\t", index=False)

    by_comparison_records = []
    for _, row in primary_df.iterrows():
        rec = {
            "tool": row["tool"],
            "modification_type": "all",
            "comparison": row["comparison"],
            "window": 0,
            "precision": row["precision"],
            "recall": row["recall"],
            "f1": row["f1"],
            "tp": 0,
            "fp": 0,
            "fn": 0,
            "tn": 0,
            "specificity": np.nan,
            "mcc": np.nan,
            "auprc": row["prauc"],
            "auroc": row["auroc"],
            "called_sites": np.nan,
            "total_truth": total_truth,
            "total_negative": np.nan,
        }
        by_comparison_records.append(rec)

    by_comparison_df = pd.DataFrame(by_comparison_records, columns=ACCURACY_BY_COMPARISON_COLS)
    by_comparison_df.to_csv(out_by_comparison, sep="\t", index=False)

    pd.DataFrame(columns=BY_NEGATIVE_TYPE_COLS).to_csv(out_by_negative_type, sep="\t", index=False)

    by_tool_records = []
    for tool, grp in primary_df.groupby("tool"):
        score_col = grp["score_column"].mode()
        best_score = score_col.iloc[0] if len(score_col) > 0 else np.nan
        rec = {
            "tool": tool,
            "best_score": best_score,
            "mean_auroc": grp["auroc"].mean(),
            "std_auroc": grp["auroc"].std(),
            "mean_prauc": grp["prauc"].mean(),
            "std_prauc": grp["prauc"].std(),
            "mean_f1": grp["f1"].mean(),
            "std_f1": grp["f1"].std(),
            "n_comparisons": len(grp),
        }
        by_tool_records.append(rec)

    by_tool_df = pd.DataFrame(by_tool_records, columns=BY_TOOL_COLS)
    by_tool_df.to_csv(out_by_tool, sep="\t", index=False)

    best_scores_records = []
    for tool, grp in primary_df.groupby("tool"):
        score_counts = grp["score_column"].value_counts()
        if len(score_counts) > 0:
            best_score = score_counts.index[0]
            n_best = score_counts.iloc[0]
        else:
            best_score = np.nan
            n_best = 0
        best_scores_records.append({
            "tool": tool,
            "best_score": best_score,
            "n_comparisons_best": n_best,
        })

    best_scores_df = pd.DataFrame(best_scores_records, columns=BEST_SCORES_COLS)
    best_scores_df.to_csv(out_best_scores, sep="\t", index=False)

called_sites_by_comp_records = []
for path in called_sites_files:
    comparison = comparison_from_called_sites_path(path)
    if comparison is None:
        continue
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        continue
    if df.empty:
        continue
    for _, row in df.iterrows():
        called_sites_by_comp_records.append({
            "tool": row["tool"],
            "comparison": comparison,
            "n_covered": row["n_covered"],
        })

called_sites_by_comp_df = pd.DataFrame(
    called_sites_by_comp_records,
    columns=["tool", "comparison", "n_covered"],
)
called_sites_by_comp_df.to_csv(out_called_sites_comp, sep="\t", index=False)

if not called_sites_by_comp_df.empty:
    called_sites_sum_df = (
        called_sites_by_comp_df
        .groupby("tool")
        .agg(total_covered=("n_covered", "sum"), n_comparisons=("comparison", "count"))
        .reset_index()
    )
else:
    called_sites_sum_df = pd.DataFrame(columns=["tool", "total_covered", "n_comparisons"])

called_sites_sum_df.to_csv(out_called_sites_sum, sep="\t", index=False)
