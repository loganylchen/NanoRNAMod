"""
Aggregate per-tool benchmark results into summary tables.

Logic:
1. Load best_metrics.tsv from each tool × comparison (native + fair modes)
2. For each tool, select the best score column by highest mean AUROC across comparisons
3. Using the selected score column, compute cross-comparison aggregated metrics
4. Merge coverage data from called_sites files
5. Produce multiple output views for downstream visualization

Inputs (from Snakemake rule benchmark_aggregate):
    native_metrics: list of best_metrics.tsv paths (native evaluation)
    fair_metrics: list of best_metrics.tsv paths (fair evaluation)
    called_sites: list of called_sites.tsv paths (per-comparison)
    truth_set: path to truth set TSV

Outputs:
    summary: Full accuracy table (all tools × comparisons, using selected score column)
    overall: Cross-comparison aggregated metrics per tool
    by_comparison: Per-comparison metrics per tool
    by_negative_type: Metrics by negative site type (from truth set labels)
    by_tool: Final per-tool summary with best score, metrics, and coverage
    best_scores: Best score column selection rationale per tool
    called_sites_comp: Coverage by tool × comparison (for colored bar chart)
    called_sites_sum: Coverage summary per tool
"""

import os
import pandas as pd
import numpy as np


def to_list(x):
    if isinstance(x, str):
        return [x]
    return list(x)


def tool_comparison_from_path(path, mode):
    """Extract tool and comparison from path like .../per_tool/{tool}/{comp}/{mode}/best_metrics.tsv"""
    parts = path.replace("\\", "/").split("/")
    # New layout: per_tool/{tool}/{comp_or_sample}/{mode}/file.tsv
    for i, part in enumerate(parts):
        if part == "per_tool" and i + 3 < len(parts):
            tool = parts[i + 1]
            comparison = parts[i + 2]
            found_mode = parts[i + 3]
            if found_mode == mode:
                return tool, comparison
    # Legacy fallback: .../native/{tool}/{comp}/file.tsv
    for i, part in enumerate(parts):
        if part == mode and i + 2 < len(parts):
            return parts[i + 1], parts[i + 2]
    return None, None


def comparison_from_called_sites_path(path):
    """Extract comparison from path like .../coverage/{comparison}/called_sites.tsv"""
    parts = path.replace("\\", "/").split("/")
    for i, part in enumerate(parts):
        if part == "coverage" and i + 1 < len(parts):
            return parts[i + 1]
    return None


def load_all_metrics(files, mode):
    """
    Load ALL score_comparison.tsv rows (not just best_metrics) for full score column evaluation.
    Falls back to best_metrics.tsv (single best row per comparison).
    """
    records = []
    for path in files:
        tool, comparison = tool_comparison_from_path(path, mode)
        if tool is None:
            continue
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            continue
        if df.empty:
            continue
        # best_metrics.tsv has one row (the best score column for that comparison)
        # We load all rows if available
        for _, row in df.iterrows():
            records.append({
                "tool": tool,
                "comparison": comparison,
                "mode": mode,
                "score_column": row.get("score_column", np.nan),
                "is_pvalue": row.get("is_pvalue", False),
                "transform": row.get("transform", "none"),
                "auroc": pd.to_numeric(row.get("auroc", np.nan), errors="coerce"),
                "prauc": pd.to_numeric(row.get("prauc", np.nan), errors="coerce"),
                "best_threshold": pd.to_numeric(row.get("best_threshold", np.nan), errors="coerce"),
                "best_threshold_original": pd.to_numeric(
                    row.get("best_threshold_original", np.nan), errors="coerce"
                ),
                "f1": pd.to_numeric(row.get("f1", np.nan), errors="coerce"),
                "precision": pd.to_numeric(row.get("precision", np.nan), errors="coerce"),
                "recall": pd.to_numeric(row.get("recall", np.nan), errors="coerce"),
                "youden_threshold": pd.to_numeric(row.get("youden_threshold", np.nan), errors="coerce"),
                "youden_threshold_original": pd.to_numeric(
                    row.get("youden_threshold_original", np.nan), errors="coerce"
                ),
                "youden_j": pd.to_numeric(row.get("youden_j", np.nan), errors="coerce"),
                "youden_sensitivity": pd.to_numeric(row.get("youden_sensitivity", np.nan), errors="coerce"),
                "youden_specificity": pd.to_numeric(row.get("youden_specificity", np.nan), errors="coerce"),
                "n_sites": pd.to_numeric(row.get("n_sites", np.nan), errors="coerce"),
                "n_positive": pd.to_numeric(row.get("n_positive", np.nan), errors="coerce"),
                "n_negative": pd.to_numeric(row.get("n_negative", np.nan), errors="coerce"),
            })
    return pd.DataFrame(records)


def select_best_score_per_tool(metrics_df):
    """
    For each tool, select the score column with the highest mean AUROC across comparisons.

    Returns:
        - selected: DataFrame filtered to only the best score column per tool
        - rationale: DataFrame explaining the selection (all candidates and their mean AUROCs)
    """
    if metrics_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    # Compute mean AUROC per tool × score_column across comparisons/samples
    candidate_stats = (
        metrics_df
        .groupby(["tool", "score_column"])
        .agg(
            mean_auroc=("auroc", "mean"),
            std_auroc=("auroc", "std"),
            mean_prauc=("prauc", "mean"),
            mean_f1=("f1", "mean"),
            mean_precision=("precision", "mean"),
            mean_recall=("recall", "mean"),
            n_comparisons=("comparison", "nunique"),
        )
        .reset_index()
    )
    # Fill NaN std with 0 (single comparison/sample tools)
    candidate_stats["std_auroc"] = candidate_stats["std_auroc"].fillna(0.0)

    # Selection criterion: mean_auroc - std_auroc (penalizes inconsistency)
    # A column with 0.75±0.01 beats 0.76±0.15
    candidate_stats["selection_score"] = (
        candidate_stats["mean_auroc"] - candidate_stats["std_auroc"]
    )

    # For each tool, pick the score column with highest selection score
    best_idx = candidate_stats.groupby("tool")["selection_score"].idxmax()
    best_selections = candidate_stats.loc[best_idx][["tool", "score_column"]].copy()
    best_selections = best_selections.rename(columns={"score_column": "best_score_column"})

    # Filter original metrics to only the best score column per tool
    selected = metrics_df.merge(
        best_selections,
        left_on=["tool", "score_column"],
        right_on=["tool", "best_score_column"],
    ).drop(columns=["best_score_column"])

    # Rationale: mark which was selected
    rationale = candidate_stats.merge(
        best_selections, on="tool", how="left"
    )
    rationale["selected"] = rationale["score_column"] == rationale["best_score_column"]
    rationale = rationale.sort_values(
        ["tool", "mean_auroc"], ascending=[True, False]
    ).reset_index(drop=True)

    return selected, rationale


def load_called_sites(files):
    """Load called_sites.tsv files, annotate with comparison."""
    records = []
    for path in files:
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
            records.append({
                "tool": row["tool"],
                "comparison": comparison,
                "n_covered": row["n_covered"],
            })
    return pd.DataFrame(records, columns=["tool", "comparison", "n_covered"])


# ===========================================================================
# Main
# ===========================================================================

# Full score_comparison files have ALL score columns evaluated per comparison
# best_metrics files have only the single best per comparison (used as fallback)
native_all_files = to_list(snakemake.input.native_all_scores)
fair_all_files = to_list(snakemake.input.fair_all_scores)
native_best_files = to_list(snakemake.input.native_metrics)
fair_best_files = to_list(snakemake.input.fair_metrics)
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

os.makedirs(os.path.dirname(out_summary), exist_ok=True)

# Load truth set
truth = pd.read_csv(truth_set_path, sep="\t")
if "label" in truth.columns:
    truth_pos = truth[truth["label"] != "-"]
else:
    truth_pos = truth
total_truth = len(truth_pos)

# Load full score evaluation (all score columns per comparison)
# Fall back to best_metrics if score_comparison files fail
native_df = load_all_metrics(native_all_files, "native")
if native_df.empty:
    native_df = load_all_metrics(native_best_files, "native")
fair_df = load_all_metrics(fair_all_files, "fair")
if fair_df.empty:
    fair_df = load_all_metrics(fair_best_files, "fair")

# Use native as primary, fair as secondary
if not native_df.empty:
    primary_df = native_df.copy()
elif not fair_df.empty:
    primary_df = fair_df.copy()
else:
    primary_df = pd.DataFrame()

# Load coverage data (independent of metrics)
called_sites_comp_df = load_called_sites(called_sites_files)
called_sites_comp_df.to_csv(out_called_sites_comp, sep="\t", index=False)

if not called_sites_comp_df.empty:
    called_sites_sum_df = (
        called_sites_comp_df
        .groupby("tool")
        .agg(total_covered=("n_covered", "sum"), n_comparisons=("comparison", "count"))
        .reset_index()
    )
else:
    called_sites_sum_df = pd.DataFrame(columns=["tool", "total_covered", "n_comparisons"])
called_sites_sum_df.to_csv(out_called_sites_sum, sep="\t", index=False)

# --- If no metrics data, write empty outputs and exit ---
if primary_df.empty:
    for path in [out_summary, out_overall, out_by_comparison,
                 out_by_negative_type, out_by_tool, out_best_scores]:
        pd.DataFrame().to_csv(path, sep="\t", index=False)
else:
    # Step 1: Select best score column per tool (by mean AUROC)
    selected_df, rationale_df = select_best_score_per_tool(primary_df)

    # --- Output: best_scores.tsv (selection rationale) ---
    rationale_df.to_csv(out_best_scores, sep="\t", index=False)

    # --- Output: accuracy_summary.tsv (all tool × comparison rows, best score only) ---
    summary_cols = [
        "tool", "comparison", "score_column", "is_pvalue", "transform",
        "auroc", "prauc", "best_threshold", "best_threshold_original",
        "f1", "precision", "recall",
        "youden_threshold", "youden_threshold_original", "youden_j",
        "youden_sensitivity", "youden_specificity",
        "n_sites", "n_positive", "n_negative",
    ]
    available_summary_cols = [c for c in summary_cols if c in selected_df.columns]
    selected_df[available_summary_cols].to_csv(out_summary, sep="\t", index=False)

    # --- Output: accuracy_summary_by_comparison.tsv ---
    # Includes both native and fair rows (mode column) for downstream visualization
    # Filter fair data to the same best score columns selected from native
    if not fair_df.empty:
        best_selections = selected_df[["tool", "score_column"]].drop_duplicates()
        fair_selected = fair_df.merge(
            best_selections, on=["tool", "score_column"], how="inner"
        )
        by_comp_combined = pd.concat([selected_df, fair_selected], ignore_index=True)
    else:
        by_comp_combined = selected_df.copy()
    by_comp_cols = available_summary_cols + (["mode"] if "mode" in by_comp_combined.columns else [])
    by_comp_available = [c for c in by_comp_cols if c in by_comp_combined.columns]
    by_comp_combined[by_comp_available].to_csv(out_by_comparison, sep="\t", index=False)

    # --- Output: accuracy_summary_overall.tsv (cross-comparison aggregate per tool) ---
    overall_records = []
    for tool, grp in selected_df.groupby("tool"):
        score_col = grp["score_column"].iloc[0]
        rec = {
            "tool": tool,
            "score_column": score_col,
            "is_pvalue": grp["is_pvalue"].iloc[0],
            "n_comparisons": len(grp),
            "auroc_mean": grp["auroc"].mean(),
            "auroc_std": grp["auroc"].std(),
            "prauc_mean": grp["prauc"].mean(),
            "prauc_std": grp["prauc"].std(),
            "f1_mean": grp["f1"].mean(),
            "f1_std": grp["f1"].std(),
            "precision_mean": grp["precision"].mean(),
            "precision_std": grp["precision"].std(),
            "recall_mean": grp["recall"].mean(),
            "recall_std": grp["recall"].std(),
            "threshold_mean": grp["best_threshold"].mean(),
            "threshold_std": grp["best_threshold"].std(),
            "youden_j_mean": grp["youden_j"].mean() if "youden_j" in grp.columns else np.nan,
            "youden_threshold_mean": grp["youden_threshold"].mean() if "youden_threshold" in grp.columns else np.nan,
            "n_sites_mean": grp["n_sites"].mean() if "n_sites" in grp.columns else np.nan,
            "n_positive_mean": grp["n_positive"].mean() if "n_positive" in grp.columns else np.nan,
            "n_negative_mean": grp["n_negative"].mean() if "n_negative" in grp.columns else np.nan,
            "total_truth": total_truth,
        }
        overall_records.append(rec)

    overall_df = pd.DataFrame(overall_records)
    overall_df = overall_df.sort_values("auroc_mean", ascending=False).reset_index(drop=True)
    overall_df.to_csv(out_overall, sep="\t", index=False)

    # --- Output: by_tool.tsv (final per-tool summary with coverage) ---
    # This is the KEY output for cross-tool comparison
    by_tool_records = []
    for tool, grp in selected_df.groupby("tool"):
        rec = {
            "tool": tool,
            "best_score_column": grp["score_column"].iloc[0],
            "is_pvalue": grp["is_pvalue"].iloc[0],
            "transform": grp["transform"].iloc[0],
            "n_comparisons": len(grp),
            "auroc_mean": round(grp["auroc"].mean(), 4),
            "auroc_std": round(grp["auroc"].std(), 4),
            "prauc_mean": round(grp["prauc"].mean(), 4),
            "f1_mean": round(grp["f1"].mean(), 4),
            "precision_mean": round(grp["precision"].mean(), 4),
            "recall_mean": round(grp["recall"].mean(), 4),
            "threshold_mean": round(grp["best_threshold"].mean(), 4),
            "threshold_original_mean": round(grp["best_threshold_original"].mean(), 6),
            "youden_j_mean": round(grp["youden_j"].mean(), 4) if "youden_j" in grp.columns and grp["youden_j"].notna().any() else np.nan,
            "youden_threshold_mean": round(grp["youden_threshold"].mean(), 4) if "youden_threshold" in grp.columns and grp["youden_threshold"].notna().any() else np.nan,
            "n_sites_mean": int(round(grp["n_sites"].mean())) if "n_sites" in grp.columns and grp["n_sites"].notna().any() else np.nan,
        }

        # Add per-comparison AUROC breakdown
        for _, comp_row in grp.iterrows():
            comp_key = f"auroc_{comp_row['comparison']}"
            rec[comp_key] = round(comp_row["auroc"], 4)

        # Merge coverage data
        tool_coverage = called_sites_comp_df[called_sites_comp_df["tool"] == tool]
        if not tool_coverage.empty:
            rec["total_sites_covered"] = int(tool_coverage["n_covered"].sum())
            rec["mean_sites_covered"] = round(tool_coverage["n_covered"].mean(), 1)
            # Per-comparison coverage
            for _, cov_row in tool_coverage.iterrows():
                cov_key = f"covered_{cov_row['comparison']}"
                rec[cov_key] = int(cov_row["n_covered"])
        else:
            rec["total_sites_covered"] = 0
            rec["mean_sites_covered"] = 0.0

        by_tool_records.append(rec)

    by_tool_df = pd.DataFrame(by_tool_records)
    by_tool_df = by_tool_df.sort_values("auroc_mean", ascending=False).reset_index(drop=True)
    by_tool_df.to_csv(out_by_tool, sep="\t", index=False)

    # --- Output: accuracy_summary_by_negative_type.tsv ---
    # Only meaningful if truth set has modification type labels
    if "label" in truth.columns:
        label_counts = truth["label"].value_counts().reset_index()
        label_counts.columns = ["label", "n_sites"]
        label_counts.to_csv(out_by_negative_type, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["label", "n_sites"]).to_csv(
            out_by_negative_type, sep="\t", index=False
        )
