"""
Generate per-tool × comparison benchmark figures.

For each tool and comparison, produces 4 figures with co-located source data:
  - roc_curve.pdf — native + fair ROC overlaid, with site counts in legend
  - pr_curve.pdf — native + fair PR overlaid, with site counts in legend
  - score_distribution.pdf — score histogram by label (pos/neg) for native mode
  - native_vs_fair.pdf — bar chart comparing metrics between modes + site count annotation

Inputs (via Snakemake):
    native_scores: score_comparison.tsv from native mode
    fair_scores: score_comparison.tsv from fair mode
    results: raw tool results TSV
    union_predictions: union of all tool-reported sites (for fair mode)
    truth_set: truth set TSV

Outputs:
    8 files in per_tool/{tool}/{comp}/figures/ (4 PDFs + 4 data TSVs)
"""

import os
import pandas as pd
import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("matplotlib is required for per-tool figures")

from sklearn.metrics import roc_curve, precision_recall_curve, auc


def load_scores(path):
    """Load score_comparison.tsv, return DataFrame or empty."""
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return pd.DataFrame()
        return df
    except Exception:
        return pd.DataFrame()


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


INF_CAP = 1000.0


def build_labels_and_scores(preds_df, truth_pos, score_col, is_pval, window, eval_sites=None, fill_value=None):
    """Build aligned labels and scores arrays.

    Args:
        preds_df: Predictions DataFrame (normalized columns: transcript, position)
        truth_pos: Truth positives DataFrame
        score_col: Name of the score column
        is_pval: Whether scores are p-values
        window: Positional window for truth matching
        eval_sites: If provided (fair mode), evaluate on these (tx, pos) pairs
        fill_value: Value for sites not in preds_df (fair mode only)

    Returns:
        (labels, scores, n_sites) or (None, None, 0) on failure
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
        # Determine label
        label = 0
        if tx in truth_dict:
            for tp in truth_dict[tx]:
                if abs(pos - tp) <= window:
                    label = 1
                    break

        # Get score
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

    # Transform p-values
    if is_pval:
        scores = -np.log10(np.clip(scores, 1e-300, None))

    # Cap inf
    scores = np.where(np.isposinf(scores), INF_CAP, scores)
    scores = np.where(np.isneginf(scores), -INF_CAP, scores)

    # Auto-detect inverted score direction: if AUROC < 0.5, negate scores
    try:
        from sklearn.metrics import roc_auc_score
        auroc_check = roc_auc_score(labels, scores)
        if auroc_check < 0.5:
            scores = -scores
    except Exception:
        pass

    return labels, scores, len(labels)


def compute_curves(labels, scores):
    """Compute ROC and PR curve data.

    Returns:
        (roc_df, roc_auc), (pr_df, pr_auc) or None on failure
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


# ===========================================================================
# Main
# ===========================================================================

native_scores_path = snakemake.input.native_scores
fair_scores_path = snakemake.input.fair_scores
results_path = snakemake.input.results
union_predictions_path = snakemake.input.union_predictions
truth_path = snakemake.input.truth_set

window_param = snakemake.params.window
if isinstance(window_param, list):
    window = max(int(w) for w in window_param)
else:
    window = int(window_param)

out_dir = os.path.dirname(snakemake.output.roc_pdf)
os.makedirs(out_dir, exist_ok=True)

native_df = load_scores(native_scores_path)
fair_df = load_scores(fair_scores_path)

# Determine best score column from native mode
best_score_col = None
is_pval = False
if not native_df.empty:
    best_row = native_df.iloc[0]
    best_score_col = best_row.get("score_column", None)
    is_pval = bool(best_row.get("is_pvalue", False))
elif not fair_df.empty:
    best_row = fair_df.iloc[0]
    best_score_col = best_row.get("score_column", None)
    is_pval = bool(best_row.get("is_pvalue", False))

# Load raw data
sep = "," if results_path.endswith(".csv") else "\t"
try:
    preds = pd.read_csv(results_path, sep=sep)
    preds = normalize_columns(preds)
    preds["position"] = pd.to_numeric(preds["position"], errors="coerce")
    preds = preds.dropna(subset=["position"])
    preds["position"] = preds["position"].astype(int)
except Exception:
    preds = pd.DataFrame()

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
try:
    union_df = pd.read_csv(union_predictions_path, sep="\t")
    union_df = normalize_columns(union_df)
    union_df["position"] = pd.to_numeric(union_df["position"], errors="coerce")
    union_df = union_df.dropna(subset=["position"])
    union_df["position"] = union_df["position"].astype(int)
    fair_eval_sites = list(zip(union_df["transcript"], union_df["position"]))
except Exception:
    fair_eval_sites = []

# Compute fair mode fill value
fair_fill = None
if is_pval:
    fair_fill = 1.0
elif best_score_col and not preds.empty and best_score_col in preds.columns:
    vals = pd.to_numeric(preds[best_score_col], errors="coerce").dropna()
    if len(vals) > 0:
        fair_fill = float(vals.min()) - 1.0
    else:
        fair_fill = -1.0

# Build native and fair labels/scores
native_labels, native_scores, n_native = None, None, 0
fair_labels, fair_scores_arr, n_fair = None, None, 0

if best_score_col and not preds.empty:
    native_labels, native_scores, n_native = build_labels_and_scores(
        preds, truth_pos, best_score_col, is_pval, window
    )
    if fair_eval_sites:
        fair_labels, fair_scores_arr, n_fair = build_labels_and_scores(
            preds, truth_pos, best_score_col, is_pval, window,
            eval_sites=fair_eval_sites, fill_value=fair_fill
        )

# Compute curves
native_roc, native_pr = compute_curves(native_labels, native_scores)
fair_roc, fair_pr = compute_curves(fair_labels, fair_scores_arr)

# Count positives/negatives for annotation
n_native_pos = int(np.sum(native_labels == 1)) if native_labels is not None else 0
n_native_neg = int(np.sum(native_labels == 0)) if native_labels is not None else 0
n_fair_pos = int(np.sum(fair_labels == 1)) if fair_labels is not None else 0
n_fair_neg = int(np.sum(fair_labels == 0)) if fair_labels is not None else 0

# --- Fig 1: ROC curves (native + fair overlaid) ---
fig, ax = plt.subplots(figsize=(5, 5))
roc_data_rows = []

if native_roc:
    roc_df, roc_auc_val = native_roc
    label = f"Native (AUC={roc_auc_val:.3f}, n={n_native})"
    ax.plot(roc_df["fpr"], roc_df["tpr"], label=label, color="#0072B2", linewidth=1.5)
    roc_df_out = roc_df.copy()
    roc_df_out["mode"] = "native"
    roc_df_out["auroc"] = roc_auc_val
    roc_df_out["n_sites"] = n_native
    roc_data_rows.append(roc_df_out)

if fair_roc:
    roc_df, roc_auc_val = fair_roc
    label = f"Fair (AUC={roc_auc_val:.3f}, n={n_fair})"
    ax.plot(roc_df["fpr"], roc_df["tpr"], label=label, color="#D55E00", linewidth=1.5, linestyle="--")
    roc_df_out = roc_df.copy()
    roc_df_out["mode"] = "fair"
    roc_df_out["auroc"] = roc_auc_val
    roc_df_out["n_sites"] = n_fair
    roc_data_rows.append(roc_df_out)

ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("ROC Curve")
ax.legend(loc="lower right", fontsize=8)
fig.tight_layout()
fig.savefig(snakemake.output.roc_pdf, dpi=300)
plt.close(fig)

if roc_data_rows:
    pd.concat(roc_data_rows).to_csv(snakemake.output.roc_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["fpr", "tpr", "mode", "auroc", "n_sites"]).to_csv(
        snakemake.output.roc_data, sep="\t", index=False)

# --- Fig 2: PR curves (native + fair overlaid) ---
fig, ax = plt.subplots(figsize=(5, 5))
pr_data_rows = []

if native_pr:
    pr_df, pr_auc_val = native_pr
    label = f"Native (AUC={pr_auc_val:.3f}, n={n_native})"
    ax.plot(pr_df["recall"], pr_df["precision"], label=label, color="#0072B2", linewidth=1.5)
    pr_df_out = pr_df.copy()
    pr_df_out["mode"] = "native"
    pr_df_out["prauc"] = pr_auc_val
    pr_df_out["n_sites"] = n_native
    pr_data_rows.append(pr_df_out)

if fair_pr:
    pr_df, pr_auc_val = fair_pr
    label = f"Fair (AUC={pr_auc_val:.3f}, n={n_fair})"
    ax.plot(pr_df["recall"], pr_df["precision"], label=label, color="#D55E00", linewidth=1.5, linestyle="--")
    pr_df_out = pr_df.copy()
    pr_df_out["mode"] = "fair"
    pr_df_out["prauc"] = pr_auc_val
    pr_df_out["n_sites"] = n_fair
    pr_data_rows.append(pr_df_out)

ax.set_xlabel("Recall")
ax.set_ylabel("Precision")
ax.set_title("Precision-Recall Curve")
ax.legend(loc="lower left", fontsize=8)
fig.tight_layout()
fig.savefig(snakemake.output.pr_pdf, dpi=300)
plt.close(fig)

if pr_data_rows:
    pd.concat(pr_data_rows).to_csv(snakemake.output.pr_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["recall", "precision", "mode", "prauc", "n_sites"]).to_csv(
        snakemake.output.pr_data, sep="\t", index=False)

# --- Fig 3: Score distribution by label ---
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
dist_data = pd.DataFrame()

for ax_idx, (mode_name, mode_labels, mode_scores, n_sites) in enumerate([
    ("Native", native_labels, native_scores, n_native),
    ("Fair", fair_labels, fair_scores_arr, n_fair),
]):
    ax = axes[ax_idx]
    if mode_labels is not None and mode_scores is not None:
        mode_dist = pd.DataFrame({"score": mode_scores, "label": mode_labels, "mode": mode_name})
        dist_data = pd.concat([dist_data, mode_dist], ignore_index=True)

        pos_scores = mode_scores[mode_labels == 1]
        neg_scores = mode_scores[mode_labels == 0]
        if len(pos_scores) > 0:
            ax.hist(pos_scores, bins=30, alpha=0.6, label=f"Positive (n={len(pos_scores)})", color="#D55E00")
        if len(neg_scores) > 0:
            ax.hist(neg_scores, bins=30, alpha=0.6, label=f"Negative (n={len(neg_scores)})", color="#0072B2")
        ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)

    score_label = f"-log10({best_score_col})" if is_pval else (best_score_col or "Score")
    ax.set_xlabel(score_label)
    ax.set_title(f"{mode_name} (n={n_sites})")

axes[0].set_ylabel("Count")
fig.suptitle("Score Distribution by Label", fontsize=12)
fig.tight_layout()
fig.savefig(snakemake.output.dist_pdf, dpi=300)
plt.close(fig)

if not dist_data.empty:
    dist_data.to_csv(snakemake.output.dist_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["score", "label", "mode"]).to_csv(snakemake.output.dist_data, sep="\t", index=False)

# --- Fig 4: Native vs Fair metric comparison + site counts ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [3, 1]})
comparison_data = []

metrics = ["auroc", "prauc", "f1", "precision", "recall"]
for mode_name, mode_df in [("native", native_df), ("fair", fair_df)]:
    if mode_df.empty:
        continue
    best = mode_df.iloc[0]
    for m in metrics:
        val = best.get(m, np.nan)
        try:
            val = float(val)
        except (ValueError, TypeError):
            val = np.nan
        comparison_data.append({"mode": mode_name, "metric": m, "value": val})

# Add site count info to the comparison data
site_info = pd.DataFrame([
    {"mode": "native", "n_sites": n_native, "n_positive": n_native_pos, "n_negative": n_native_neg},
    {"mode": "fair", "n_sites": n_fair, "n_positive": n_fair_pos, "n_negative": n_fair_neg},
])

if comparison_data:
    comp_df = pd.DataFrame(comparison_data)
    # Merge site info into output
    full_output = comp_df.merge(site_info, on="mode", how="left")
    full_output.to_csv(snakemake.output.nvf_data, sep="\t", index=False)

    native_vals = comp_df[comp_df["mode"] == "native"].set_index("metric")["value"]
    fair_vals = comp_df[comp_df["mode"] == "fair"].set_index("metric")["value"]

    x = np.arange(len(metrics))
    width = 0.35
    if not native_vals.empty:
        ax1.bar(x - width / 2, [native_vals.get(m, 0) for m in metrics], width,
                label=f"Native (n={n_native})", color="#0072B2")
    if not fair_vals.empty:
        ax1.bar(x + width / 2, [fair_vals.get(m, 0) for m in metrics], width,
                label=f"Fair (n={n_fair})", color="#D55E00")

    ax1.set_xticks(x)
    ax1.set_xticklabels([m.upper() for m in metrics], rotation=45, ha="right")
    ax1.set_ylabel("Score")
    ax1.set_title("Native vs Fair Mode Comparison")
    ax1.legend(fontsize=8)
    ax1.set_ylim(0, 1.1)
else:
    pd.DataFrame(columns=["mode", "metric", "value", "n_sites", "n_positive", "n_negative"]).to_csv(
        snakemake.output.nvf_data, sep="\t", index=False)
    ax1.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax1.transAxes)

# Right panel: site count breakdown
site_labels = ["Native", "Fair"]
site_pos = [n_native_pos, n_fair_pos]
site_neg = [n_native_neg, n_fair_neg]
x_bar = np.arange(2)
bar_w = 0.35

ax2.bar(x_bar - bar_w / 2, site_pos, bar_w, label="Positive", color="#D55E00", alpha=0.8)
ax2.bar(x_bar + bar_w / 2, site_neg, bar_w, label="Negative", color="#0072B2", alpha=0.8)
ax2.set_xticks(x_bar)
ax2.set_xticklabels(site_labels)
ax2.set_ylabel("Count")
ax2.set_title("Evaluation Sites")
ax2.legend(fontsize=7)

# Annotate total counts
for i, (total, pos, neg) in enumerate(zip([n_native, n_fair], site_pos, site_neg)):
    ax2.text(i, max(pos, neg) * 1.05, f"n={total}", ha="center", fontsize=8, fontweight="bold")

fig.tight_layout()
fig.savefig(snakemake.output.nvf_pdf, dpi=300)
plt.close(fig)
