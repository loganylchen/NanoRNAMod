"""
Generate cross-tool comparison figures: overlaid ROC/PR curves, AUROC/AUPRC
bar chart, and score distribution violin/box plots.

For each comparison, all active tools are overlaid on the same plot so their
discriminative performance can be compared directly.

Inputs (via Snakemake):
    results: list of all tool result TSVs
    truth_set: truth set TSV
    best_scores: cross-tool best_scores.tsv (best score column per tool)
    union_predictions: per-comparison union prediction files

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


def build_labels_and_scores(preds_df, truth_pos, score_col, is_pval, window):
    """Build aligned labels and scores arrays (native mode).

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

    # Evaluate on sites the tool reports that are on truth transcripts
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
                continue
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

# Process each tool
tool_curves = {}  # tool -> {roc: (df, auc), pr: (df, auc), labels, scores}

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

    labels, scores, n_sites = build_labels_and_scores(
        preds, truth_pos, score_col, is_pval, window
    )
    if labels is None:
        continue

    roc_result, pr_result = compute_curves(labels, scores)

    # For tools with multiple comparisons, keep the one with more sites
    if tool in tool_curves and tool_curves[tool]["n_sites"] >= n_sites:
        continue

    tool_curves[tool] = {
        "roc": roc_result,
        "pr": pr_result,
        "labels": labels,
        "scores": scores,
        "n_sites": n_sites,
        "score_col": score_col,
        "is_pval": is_pval,
    }

# Sort tools alphabetically for consistent ordering
sorted_tools = sorted(tool_curves.keys())

# --- Fig 1: All-tools overlaid ROC curves ---
fig, ax = plt.subplots(figsize=(6, 6))
roc_data_rows = []

for idx, tool in enumerate(sorted_tools):
    info = tool_curves[tool]
    if info["roc"] is None:
        continue
    roc_df, roc_auc_val = info["roc"]
    color = tool_color(idx)
    label = f"{tool} (AUROC={roc_auc_val:.3f}, n={info['n_sites']})"
    ax.plot(roc_df["fpr"], roc_df["tpr"], label=label, color=color, linewidth=1.5)

    out_df = roc_df.copy()
    out_df["tool"] = tool
    out_df["auroc"] = roc_auc_val
    out_df["n_sites"] = info["n_sites"]
    roc_data_rows.append(out_df)

ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=0.8)
ax.set_xlabel("False Positive Rate", fontsize=8)
ax.set_ylabel("True Positive Rate", fontsize=8)
ax.set_title("Cross-Tool ROC Comparison", fontsize=10)
ax.legend(loc="lower right", fontsize=7, framealpha=0.9)
ax.tick_params(labelsize=8)
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
fig.tight_layout()
fig.savefig(snakemake.output.roc_pdf, dpi=300)
plt.close(fig)

if roc_data_rows:
    pd.concat(roc_data_rows).to_csv(snakemake.output.roc_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["fpr", "tpr", "tool", "auroc", "n_sites"]).to_csv(
        snakemake.output.roc_data, sep="\t", index=False
    )

# --- Fig 2: All-tools overlaid PR curves ---
fig, ax = plt.subplots(figsize=(6, 6))
pr_data_rows = []

for idx, tool in enumerate(sorted_tools):
    info = tool_curves[tool]
    if info["pr"] is None:
        continue
    pr_df, pr_auc_val = info["pr"]
    color = tool_color(idx)
    label = f"{tool} (AUPRC={pr_auc_val:.3f}, n={info['n_sites']})"
    ax.plot(pr_df["recall"], pr_df["precision"], label=label, color=color, linewidth=1.5)

    out_df = pr_df.copy()
    out_df["tool"] = tool
    out_df["auprc"] = pr_auc_val
    out_df["n_sites"] = info["n_sites"]
    pr_data_rows.append(out_df)

# Baseline: prevalence line
all_labels = np.concatenate([tool_curves[t]["labels"] for t in sorted_tools if tool_curves[t]["labels"] is not None])
if len(all_labels) > 0:
    prevalence = np.mean(all_labels)
    ax.axhline(y=prevalence, color="k", linestyle="--", alpha=0.3, linewidth=0.8,
               label=f"Baseline ({prevalence:.3f})")

ax.set_xlabel("Recall", fontsize=8)
ax.set_ylabel("Precision", fontsize=8)
ax.set_title("Cross-Tool Precision-Recall Comparison", fontsize=10)
ax.legend(loc="lower left", fontsize=7, framealpha=0.9)
ax.tick_params(labelsize=8)
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.05)
fig.tight_layout()
fig.savefig(snakemake.output.pr_pdf, dpi=300)
plt.close(fig)

if pr_data_rows:
    pd.concat(pr_data_rows).to_csv(snakemake.output.pr_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["recall", "precision", "tool", "auprc", "n_sites"]).to_csv(
        snakemake.output.pr_data, sep="\t", index=False
    )

# --- Fig 3: AUROC/AUPRC bar chart ---
fig, ax = plt.subplots(figsize=(max(4, len(sorted_tools) * 1.2), 5))
bar_data_rows = []

auroc_vals = []
auprc_vals = []
tool_labels = []

for tool in sorted_tools:
    info = tool_curves[tool]
    auroc_val = info["roc"][1] if info["roc"] else np.nan
    auprc_val = info["pr"][1] if info["pr"] else np.nan
    auroc_vals.append(auroc_val)
    auprc_vals.append(auprc_val)
    tool_labels.append(tool)
    bar_data_rows.append({
        "tool": tool,
        "auroc": auroc_val,
        "auprc": auprc_val,
        "n_sites": info["n_sites"],
    })

x = np.arange(len(tool_labels))
width = 0.35

ax.bar(x - width / 2, auroc_vals, width, label="AUROC", color="#0072B2", alpha=0.85)
ax.bar(x + width / 2, auprc_vals, width, label="AUPRC", color="#D55E00", alpha=0.85)

# Value labels on bars
for i, (a, p) in enumerate(zip(auroc_vals, auprc_vals)):
    if not np.isnan(a):
        ax.text(i - width / 2, a + 0.01, f"{a:.2f}", ha="center", va="bottom", fontsize=6)
    if not np.isnan(p):
        ax.text(i + width / 2, p + 0.01, f"{p:.2f}", ha="center", va="bottom", fontsize=6)

ax.set_xticks(x)
ax.set_xticklabels(tool_labels, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Score", fontsize=8)
ax.set_title("AUROC / AUPRC Comparison", fontsize=10)
ax.legend(fontsize=8)
ax.set_ylim(0, 1.15)
ax.tick_params(labelsize=8)
fig.tight_layout()
fig.savefig(snakemake.output.bar_pdf, dpi=300)
plt.close(fig)

if bar_data_rows:
    pd.DataFrame(bar_data_rows).to_csv(snakemake.output.bar_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["tool", "auroc", "auprc", "n_sites"]).to_csv(
        snakemake.output.bar_data, sep="\t", index=False
    )

# --- Fig 4: Score distribution violin/box plots (TP vs FP per tool) ---
fig, ax = plt.subplots(figsize=(max(5, len(sorted_tools) * 1.5), 5))
dist_data_rows = []

positions = []
violin_data_pos = []
violin_data_neg = []
pos_idx = 1

for idx, tool in enumerate(sorted_tools):
    info = tool_curves[tool]
    labels = info["labels"]
    scores = info["scores"]
    if labels is None or scores is None:
        continue

    score_label = f"-log10({info['score_col']})" if info["is_pval"] else info["score_col"]

    pos_scores = scores[labels == 1]
    neg_scores = scores[labels == 0]

    for s, lbl in [(pos_scores, "TP"), (neg_scores, "FP")]:
        for v in s:
            dist_data_rows.append({"tool": tool, "label": lbl, "score": v, "score_column": score_label})

    # Box plots side by side
    bp_pos = pos_idx
    if len(pos_scores) > 1:
        parts = ax.boxplot(
            [pos_scores], positions=[bp_pos - 0.2], widths=0.3,
            patch_artist=True, showfliers=False,
            boxprops=dict(facecolor="#D55E00", alpha=0.6),
            medianprops=dict(color="black"),
        )
    if len(neg_scores) > 1:
        parts = ax.boxplot(
            [neg_scores], positions=[bp_pos + 0.2], widths=0.3,
            patch_artist=True, showfliers=False,
            boxprops=dict(facecolor="#0072B2", alpha=0.6),
            medianprops=dict(color="black"),
        )
    positions.append(bp_pos)
    pos_idx += 1

ax.set_xticks(positions)
ax.set_xticklabels(sorted_tools, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Score (transformed)", fontsize=8)
ax.set_title("Score Distribution: TP vs FP", fontsize=10)
ax.tick_params(labelsize=8)

# Custom legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#D55E00", alpha=0.6, label="True Positive"),
    Patch(facecolor="#0072B2", alpha=0.6, label="False Positive"),
]
ax.legend(handles=legend_elements, fontsize=8, loc="upper right")

fig.tight_layout()
fig.savefig(snakemake.output.dist_pdf, dpi=300)
plt.close(fig)

if dist_data_rows:
    pd.DataFrame(dist_data_rows).to_csv(snakemake.output.dist_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["tool", "label", "score", "score_column"]).to_csv(
        snakemake.output.dist_data, sep="\t", index=False
    )
