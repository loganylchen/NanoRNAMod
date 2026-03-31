"""
Generate cross-tool comparison figures: overlaid ROC/PR curves, AUROC/AUPRC
bar chart, and score distribution box plots.

Both native and fair evaluation modes are shown on each plot so their
discriminative performance can be compared directly.

Inputs (via Snakemake):
    results: list of all tool result TSVs
    truth_set: truth set TSV
    best_scores: cross-tool best_scores.tsv (best score column per tool)
    union_predictions: per-comparison union prediction files (for fair mode)

Outputs:
    4 PDFs + 4 source-data TSVs in cross_tool/figures/
"""

import os
import re
import pandas as pd
import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("matplotlib is required for cross-tool figures")

from sklearn.metrics import roc_curve, precision_recall_curve, auc


# ---------------------------------------------------------------------------
# Helpers (aligned with benchmark_per_tool_figures.py)
# ---------------------------------------------------------------------------

INF_CAP = 1000.0

# Okabe-Ito colorblind-friendly palette (matches R figures)
OKABE_ITO = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#999999",  # grey
]


def tool_color(idx):
    """Return a color from the Okabe-Ito palette, cycling if needed."""
    return OKABE_ITO[idx % len(OKABE_ITO)]


def normalize_columns(df):
    """Normalize column names to standard transcript/position."""
    col_mapping = {}
    for col in ["transcript", "id", "ref_id", "chrom", "transcript_id"]:
        if col in df.columns:
            col_mapping[col] = "transcript"
            break
    for col in ["position", "pos", "start", "start_loc", "transcript_pos"]:
        if col in df.columns and col not in col_mapping.values():
            col_mapping[col] = "position"
            break
    if col_mapping:
        df = df.rename(columns=col_mapping)
    return df


def build_labels_and_scores(preds_df, truth_pos, score_col, is_pval, window,
                            eval_sites=None, fill_value=None):
    """Build aligned labels and scores arrays.

    Args:
        preds_df: Predictions DataFrame (normalized columns: transcript, position)
        truth_pos: Truth positives DataFrame
        score_col: Name of the score column
        is_pval: Whether scores are p-values
        window: Positional window for truth matching
        eval_sites: If provided (fair mode), evaluate on these (tx, pos) pairs
        fill_value: Value for sites not in preds_df (fair mode only)

    Returns (labels, scores, n_sites) or (None, None, 0) on failure.
    """
    if score_col not in preds_df.columns:
        return None, None, 0

    # Build truth dict
    truth_dict = {}
    for _, row in truth_pos.iterrows():
        tx = row["transcript"]
        if tx not in truth_dict:
            truth_dict[tx] = []
        truth_dict[tx].append(int(row["position"]))

    # Build index from predictions
    pred_index = {}
    for _, row in preds_df.iterrows():
        key = (row["transcript"], int(row["position"]))
        try:
            val = float(row[score_col])
        except (ValueError, TypeError):
            continue
        pred_index[key] = val

    if eval_sites is None:
        # Native mode: evaluate on sites the tool reports that are on truth transcripts
        truth_txs = set(truth_pos["transcript"].unique())
        eval_sites = [(tx, pos) for (tx, pos) in pred_index.keys() if tx in truth_txs]

    labels = []
    scores = []
    for tx, pos in eval_sites:
        pos = int(pos)
        label = 0
        if tx in truth_dict:
            for tp in truth_dict[tx]:
                if abs(pos - tp) <= window:
                    label = 1
                    break
        key = (tx, pos)
        if key in pred_index:
            s = pred_index[key]
            if np.isnan(s):
                if fill_value is not None:
                    s = fill_value
                else:
                    continue
        elif fill_value is not None:
            s = fill_value
        else:
            continue
        labels.append(label)
        scores.append(s)

    if len(labels) < 2 or len(set(labels)) < 2:
        return None, None, len(labels)

    labels = np.array(labels)
    scores = np.array(scores)

    if is_pval:
        scores = -np.log10(np.clip(scores, 1e-300, None))

    scores = np.where(np.isposinf(scores), INF_CAP, scores)
    scores = np.where(np.isneginf(scores), -INF_CAP, scores)

    return labels, scores, len(labels)


def compute_curves(labels, scores):
    """Compute ROC and PR curve data.

    Returns (roc_df, roc_auc), (pr_df, pr_auc) or (None, None).
    """
    if labels is None or scores is None:
        return None, None

    try:
        fpr, tpr, _ = roc_curve(labels, scores)
        roc_auc_val = auc(fpr, tpr)
        roc_df = pd.DataFrame({"fpr": fpr, "tpr": tpr})
    except Exception:
        roc_df, roc_auc_val = None, None

    try:
        prec, rec, _ = precision_recall_curve(labels, scores)
        pr_auc_val = auc(rec, prec)
        pr_df = pd.DataFrame({"recall": rec, "precision": prec})
    except Exception:
        pr_df, pr_auc_val = None, None

    roc_result = (roc_df, roc_auc_val) if roc_df is not None else None
    pr_result = (pr_df, pr_auc_val) if pr_df is not None else None

    return roc_result, pr_result


def extract_tool_from_path(path):
    """Extract tool name from a results path like .../modifications/{tool}/.../{tool}_results.tsv."""
    m = re.search(r"/modifications/([^/]+)/", path)
    return m.group(1) if m else os.path.basename(path).replace("_results.tsv", "")


# ===========================================================================
# Main
# ===========================================================================

result_paths = snakemake.input.results
truth_path = snakemake.input.truth_set
best_scores_path = snakemake.input.best_scores
union_paths = snakemake.input.union_predictions
if isinstance(union_paths, str):
    union_paths = [union_paths]

window_param = snakemake.params.window
if isinstance(window_param, list):
    window = max(int(w) for w in window_param)
else:
    window = int(window_param)

out_dir = os.path.dirname(snakemake.output.roc_pdf)
os.makedirs(out_dir, exist_ok=True)

# Load best scores table
best_scores_df = pd.read_csv(best_scores_path, sep="\t")

# Build lookup: tool -> (best_score_col, is_pvalue)
tool_score_info = {}
for _, row in best_scores_df.iterrows():
    tool = row["tool"]
    score_col = row.get("best_score_column", row.get("score_column", None))
    is_pval = bool(row.get("is_pvalue", False))
    if score_col:
        tool_score_info[tool] = (score_col, is_pval)

# Load truth set
truth = pd.read_csv(truth_path, sep="\t")
truth = normalize_columns(truth)
if "label" in truth.columns:
    truth_pos = truth[truth["label"] != "-"].copy()
else:
    truth_pos = truth.copy()
truth_pos["position"] = pd.to_numeric(truth_pos["position"], errors="coerce")
truth_pos = truth_pos.dropna(subset=["position"])
truth_pos["position"] = truth_pos["position"].astype(int)

# Load union predictions for fair mode
all_fair_eval_sites = []
for upath in union_paths:
    try:
        udf = pd.read_csv(upath, sep="\t")
        udf = normalize_columns(udf)
        udf["position"] = pd.to_numeric(udf["position"], errors="coerce")
        udf = udf.dropna(subset=["position"])
        udf["position"] = udf["position"].astype(int)
        all_fair_eval_sites.extend(list(zip(udf["transcript"], udf["position"])))
    except Exception:
        continue
# De-duplicate
fair_eval_sites = list(set(all_fair_eval_sites)) if all_fair_eval_sites else []

# Load and process each tool's predictions
tool_preds = {}  # tool -> (preds_df, score_col, is_pval)

for rpath in result_paths:
    tool = extract_tool_from_path(rpath)
    if tool not in tool_score_info:
        continue

    score_col, is_pval = tool_score_info[tool]

    sep = "," if rpath.endswith(".csv") else "\t"
    try:
        preds = pd.read_csv(rpath, sep=sep)
        preds = normalize_columns(preds)
        preds["position"] = pd.to_numeric(preds["position"], errors="coerce")
        preds = preds.dropna(subset=["position"])
        preds["position"] = preds["position"].astype(int)
    except Exception:
        continue

    if preds.empty:
        continue

    # Keep the result with most rows if a tool has multiple comparisons
    if tool in tool_preds and len(tool_preds[tool][0]) >= len(preds):
        continue
    tool_preds[tool] = (preds, score_col, is_pval)

# Build curves for both modes
# mode_curves[mode][tool] = {roc, pr, labels, scores, n_sites, score_col, is_pval}
mode_curves = {"native": {}, "fair": {}}

for tool, (preds, score_col, is_pval) in tool_preds.items():
    # Native mode
    labels, scores, n_sites = build_labels_and_scores(
        preds, truth_pos, score_col, is_pval, window
    )
    if labels is not None:
        roc_result, pr_result = compute_curves(labels, scores)
        mode_curves["native"][tool] = {
            "roc": roc_result, "pr": pr_result,
            "labels": labels, "scores": scores,
            "n_sites": n_sites, "score_col": score_col, "is_pval": is_pval,
        }

    # Fair mode
    if fair_eval_sites:
        fair_fill = 1.0 if is_pval else None
        if not is_pval and score_col in preds.columns:
            vals = pd.to_numeric(preds[score_col], errors="coerce").dropna()
            if len(vals) > 0:
                fair_fill = float(vals.min())
        labels_f, scores_f, n_sites_f = build_labels_and_scores(
            preds, truth_pos, score_col, is_pval, window,
            eval_sites=fair_eval_sites, fill_value=fair_fill,
        )
        if labels_f is not None:
            roc_result_f, pr_result_f = compute_curves(labels_f, scores_f)
            mode_curves["fair"][tool] = {
                "roc": roc_result_f, "pr": pr_result_f,
                "labels": labels_f, "scores": scores_f,
                "n_sites": n_sites_f, "score_col": score_col, "is_pval": is_pval,
            }

# All tools present in either mode
all_tools = sorted(set(list(mode_curves["native"].keys()) + list(mode_curves["fair"].keys())))
has_fair = len(mode_curves["fair"]) > 0

# ---------------------------------------------------------------------------
# Fig 1: All-tools overlaid ROC curves (native solid, fair dashed)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 6))
roc_data_rows = []

for idx, tool in enumerate(all_tools):
    color = tool_color(idx)
    for mode, ls in [("native", "-"), ("fair", "--")]:
        info = mode_curves[mode].get(tool)
        if info is None or info["roc"] is None:
            continue
        roc_df, roc_auc_val = info["roc"]
        suffix = f" [{mode}]" if has_fair else ""
        label = f"{tool}{suffix} (AUROC={roc_auc_val:.3f}, n={info['n_sites']})"
        ax.plot(roc_df["fpr"], roc_df["tpr"], label=label, color=color,
                linewidth=1.5, linestyle=ls)

        out_df = roc_df.copy()
        out_df["tool"] = tool
        out_df["mode"] = mode
        out_df["auroc"] = roc_auc_val
        out_df["n_sites"] = info["n_sites"]
        roc_data_rows.append(out_df)

ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=0.8)
ax.set_xlabel("False Positive Rate", fontsize=8)
ax.set_ylabel("True Positive Rate", fontsize=8)
title_suffix = " (solid=native, dashed=fair)" if has_fair else ""
ax.set_title(f"Cross-Tool ROC Comparison{title_suffix}", fontsize=10)
ax.legend(loc="lower right", fontsize=6, framealpha=0.9)
ax.tick_params(labelsize=8)
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
fig.tight_layout()
fig.savefig(snakemake.output.roc_pdf, dpi=300)
plt.close(fig)

if roc_data_rows:
    pd.concat(roc_data_rows).to_csv(snakemake.output.roc_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["fpr", "tpr", "tool", "mode", "auroc", "n_sites"]).to_csv(
        snakemake.output.roc_data, sep="\t", index=False
    )

# ---------------------------------------------------------------------------
# Fig 2: All-tools overlaid PR curves (native solid, fair dashed)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 6))
pr_data_rows = []

for idx, tool in enumerate(all_tools):
    color = tool_color(idx)
    for mode, ls in [("native", "-"), ("fair", "--")]:
        info = mode_curves[mode].get(tool)
        if info is None or info["pr"] is None:
            continue
        pr_df, pr_auc_val = info["pr"]
        suffix = f" [{mode}]" if has_fair else ""
        label = f"{tool}{suffix} (AUPRC={pr_auc_val:.3f}, n={info['n_sites']})"
        ax.plot(pr_df["recall"], pr_df["precision"], label=label, color=color,
                linewidth=1.5, linestyle=ls)

        out_df = pr_df.copy()
        out_df["tool"] = tool
        out_df["mode"] = mode
        out_df["auprc"] = pr_auc_val
        out_df["n_sites"] = info["n_sites"]
        pr_data_rows.append(out_df)

# Baseline: prevalence line
all_native_labels = [mode_curves["native"][t]["labels"] for t in all_tools
                     if t in mode_curves["native"] and mode_curves["native"][t]["labels"] is not None]
if all_native_labels:
    concatenated = np.concatenate(all_native_labels)
    if len(concatenated) > 0:
        prevalence = np.mean(concatenated)
        ax.axhline(y=prevalence, color="k", linestyle="--", alpha=0.3, linewidth=0.8,
                   label=f"Baseline ({prevalence:.3f})")

ax.set_xlabel("Recall", fontsize=8)
ax.set_ylabel("Precision", fontsize=8)
title_suffix = " (solid=native, dashed=fair)" if has_fair else ""
ax.set_title(f"Cross-Tool Precision-Recall Comparison{title_suffix}", fontsize=10)
ax.legend(loc="lower left", fontsize=6, framealpha=0.9)
ax.tick_params(labelsize=8)
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.05)
fig.tight_layout()
fig.savefig(snakemake.output.pr_pdf, dpi=300)
plt.close(fig)

if pr_data_rows:
    pd.concat(pr_data_rows).to_csv(snakemake.output.pr_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["recall", "precision", "tool", "mode", "auprc", "n_sites"]).to_csv(
        snakemake.output.pr_data, sep="\t", index=False
    )

# ---------------------------------------------------------------------------
# Fig 3: AUROC/AUPRC grouped bar chart (native + fair side by side)
# ---------------------------------------------------------------------------
bar_data_rows = []
modes_to_plot = ["native", "fair"] if has_fair else ["native"]

for tool in all_tools:
    for mode in modes_to_plot:
        info = mode_curves[mode].get(tool)
        if info is None:
            continue
        auroc_val = info["roc"][1] if info["roc"] else np.nan
        auprc_val = info["pr"][1] if info["pr"] else np.nan
        bar_data_rows.append({
            "tool": tool, "mode": mode,
            "auroc": auroc_val, "auprc": auprc_val,
            "n_sites": info["n_sites"],
        })

if bar_data_rows:
    bar_df = pd.DataFrame(bar_data_rows)
else:
    bar_df = pd.DataFrame(columns=["tool", "mode", "auroc", "auprc", "n_sites"])

fig, ax = plt.subplots(figsize=(max(4, len(all_tools) * (1.8 if has_fair else 1.2)), 5))

if has_fair:
    # 4 bars per tool: native AUROC, native AUPRC, fair AUROC, fair AUPRC
    n = len(all_tools)
    x = np.arange(n)
    w = 0.18
    offsets = [-1.5 * w, -0.5 * w, 0.5 * w, 1.5 * w]
    labels_legend = ["AUROC (native)", "AUPRC (native)", "AUROC (fair)", "AUPRC (fair)"]
    colors = ["#0072B2", "#D55E00", "#56B4E9", "#E69F00"]

    for bar_idx, (mode, metric, col) in enumerate([
        ("native", "auroc", "#0072B2"), ("native", "auprc", "#D55E00"),
        ("fair", "auroc", "#56B4E9"), ("fair", "auprc", "#E69F00"),
    ]):
        vals = []
        for tool in all_tools:
            row = bar_df[(bar_df["tool"] == tool) & (bar_df["mode"] == mode)]
            vals.append(float(row[metric].iloc[0]) if len(row) > 0 else np.nan)
        ax.bar(x + offsets[bar_idx], vals, w, label=labels_legend[bar_idx], color=col, alpha=0.85)
        for i, v in enumerate(vals):
            if not np.isnan(v):
                ax.text(i + offsets[bar_idx], v + 0.01, f"{v:.2f}", ha="center", va="bottom", fontsize=5)
else:
    x = np.arange(len(all_tools))
    w = 0.35
    auroc_vals = []
    auprc_vals = []
    for tool in all_tools:
        row = bar_df[bar_df["tool"] == tool]
        auroc_vals.append(float(row["auroc"].iloc[0]) if len(row) > 0 else np.nan)
        auprc_vals.append(float(row["auprc"].iloc[0]) if len(row) > 0 else np.nan)
    ax.bar(x - w / 2, auroc_vals, w, label="AUROC", color="#0072B2", alpha=0.85)
    ax.bar(x + w / 2, auprc_vals, w, label="AUPRC", color="#D55E00", alpha=0.85)
    for i, (a, p) in enumerate(zip(auroc_vals, auprc_vals)):
        if not np.isnan(a):
            ax.text(i - w / 2, a + 0.01, f"{a:.2f}", ha="center", va="bottom", fontsize=6)
        if not np.isnan(p):
            ax.text(i + w / 2, p + 0.01, f"{p:.2f}", ha="center", va="bottom", fontsize=6)

ax.set_xticks(np.arange(len(all_tools)))
ax.set_xticklabels(all_tools, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Score", fontsize=8)
ax.set_title("AUROC / AUPRC Comparison", fontsize=10)
ax.legend(fontsize=7)
ax.set_ylim(0, 1.15)
ax.tick_params(labelsize=8)
fig.tight_layout()
fig.savefig(snakemake.output.bar_pdf, dpi=300)
plt.close(fig)

bar_df.to_csv(snakemake.output.bar_data, sep="\t", index=False)

# ---------------------------------------------------------------------------
# Fig 4: Score distribution box plots (TP vs FP, per-tool facets, free y-axes)
# ---------------------------------------------------------------------------
from matplotlib.patches import Patch

dist_data_rows = []

# Collect data per tool (native mode only for score distributions)
plot_tools = []
for tool in all_tools:
    info = mode_curves["native"].get(tool)
    if info is None:
        continue
    labels = info["labels"]
    scores = info["scores"]
    if labels is None or scores is None:
        continue
    score_label = f"-log10({info['score_col']})" if info["is_pval"] else info["score_col"]
    pos_scores = scores[labels == 1]
    neg_scores = scores[labels == 0]
    if len(pos_scores) < 2 and len(neg_scores) < 2:
        continue
    for s, lbl in [(pos_scores, "TP"), (neg_scores, "FP")]:
        for v in s:
            dist_data_rows.append({"tool": tool, "label": lbl, "score": v, "score_column": score_label})
    plot_tools.append((tool, score_label, pos_scores, neg_scores))

n_tools = max(len(plot_tools), 1)
ncols = min(n_tools, 4)
nrows = (n_tools + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows), squeeze=False)
axes_flat = axes.flatten()

for i, (tool, score_label, pos_scores, neg_scores) in enumerate(plot_tools):
    ax = axes_flat[i]
    box_data = []
    box_labels = []
    box_colors = []
    if len(pos_scores) > 1:
        box_data.append(pos_scores)
        box_labels.append("TP")
        box_colors.append("#D55E00")
    if len(neg_scores) > 1:
        box_data.append(neg_scores)
        box_labels.append("FP")
        box_colors.append("#0072B2")
    if box_data:
        bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True, showfliers=False,
                        medianprops=dict(color="black"))
        for patch, color in zip(bp["boxes"], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
    ax.set_title(tool, fontsize=9, fontweight="bold")
    ax.set_ylabel(score_label, fontsize=7)
    ax.tick_params(labelsize=7)

# Hide unused axes
for j in range(len(plot_tools), len(axes_flat)):
    axes_flat[j].set_visible(False)

# Shared legend
legend_elements = [
    Patch(facecolor="#D55E00", alpha=0.6, label="True Positive"),
    Patch(facecolor="#0072B2", alpha=0.6, label="False Positive"),
]
fig.legend(handles=legend_elements, fontsize=8, loc="upper right", ncol=2)
fig.suptitle("Score Distribution: TP vs FP", fontsize=11, y=1.01)

fig.tight_layout()
fig.savefig(snakemake.output.dist_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)

if dist_data_rows:
    pd.DataFrame(dist_data_rows).to_csv(snakemake.output.dist_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["tool", "label", "score", "score_column"]).to_csv(
        snakemake.output.dist_data, sep="\t", index=False
    )
