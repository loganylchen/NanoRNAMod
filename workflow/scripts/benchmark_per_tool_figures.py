"""
Generate per-tool × comparison benchmark figures.

For each tool and comparison, produces 4 figures with co-located source data:
  - roc_curve.pdf — native + fair ROC overlaid
  - pr_curve.pdf — native + fair PR overlaid
  - score_distribution.pdf — score histogram by label (pos/neg)
  - native_vs_fair.pdf — bar chart comparing metrics between modes

Inputs (via Snakemake):
    native_scores: score_comparison.tsv from native mode
    fair_scores: score_comparison.tsv from fair mode
    results: raw tool results TSV
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


def build_roc_pr_data(results_path, truth_path, score_col, is_pvalue, window):
    """Build ROC/PR curve data from raw results and truth set."""
    sep = "," if results_path.endswith(".csv") else "\t"
    try:
        preds = pd.read_csv(results_path, sep=sep)
    except Exception:
        return None, None

    preds = normalize_columns(preds)
    if "transcript" not in preds.columns or "position" not in preds.columns:
        return None, None
    if score_col not in preds.columns:
        return None, None

    truth = pd.read_csv(truth_path, sep="\t")
    truth = normalize_columns(truth)
    if "label" in truth.columns:
        truth_pos = truth[truth["label"] != "-"].copy()
    else:
        truth_pos = truth.copy()

    if truth_pos.empty:
        return None, None

    truth_pos["position"] = pd.to_numeric(truth_pos["position"], errors="coerce")
    truth_pos = truth_pos.dropna(subset=["position"])
    truth_pos["position"] = truth_pos["position"].astype(int)

    preds["position"] = pd.to_numeric(preds["position"], errors="coerce")
    preds = preds.dropna(subset=["position"])
    preds["position"] = preds["position"].astype(int)

    # Build truth dict
    truth_dict = {}
    for _, row in truth_pos.iterrows():
        tx = row["transcript"]
        if tx not in truth_dict:
            truth_dict[tx] = []
        truth_dict[tx].append(int(row["position"]))

    # Build labels for all prediction sites on truth transcripts
    truth_txs = set(truth_pos["transcript"].unique())
    pred_on_truth = preds[preds["transcript"].isin(truth_txs)].copy()

    labels = []
    scores = []
    for _, row in pred_on_truth.iterrows():
        tx = row["transcript"]
        pos = int(row["position"])
        label = 0
        if tx in truth_dict:
            for tp in truth_dict[tx]:
                if abs(pos - tp) <= window:
                    label = 1
                    break
        try:
            s = float(row[score_col])
        except (ValueError, TypeError):
            continue
        if np.isnan(s):
            continue
        labels.append(label)
        scores.append(s)

    if len(labels) < 2 or len(set(labels)) < 2:
        return None, None

    labels = np.array(labels)
    scores = np.array(scores)

    if is_pvalue:
        scores = -np.log10(np.clip(scores, 1e-300, None))

    fpr, tpr, _ = roc_curve(labels, scores)
    roc_auc = auc(fpr, tpr)

    prec, rec, _ = precision_recall_curve(labels, scores)
    pr_auc = auc(rec, prec)

    roc_data = pd.DataFrame({"fpr": fpr, "tpr": tpr})
    pr_data = pd.DataFrame({"recall": rec, "precision": prec})

    return (roc_data, roc_auc), (pr_data, pr_auc)


# ===========================================================================
# Main
# ===========================================================================

native_scores_path = snakemake.input.native_scores
fair_scores_path = snakemake.input.fair_scores
results_path = snakemake.input.results
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

# --- Fig 1: ROC curves (native + fair overlaid) ---
fig, ax = plt.subplots(figsize=(5, 5))
roc_data_rows = []

if best_score_col:
    roc_result, _ = build_roc_pr_data(results_path, truth_path, best_score_col, is_pval, window)
    if roc_result:
        roc_df, roc_auc_val = roc_result
        ax.plot(roc_df["fpr"], roc_df["tpr"], label=f"Native (AUC={roc_auc_val:.3f})", color="#0072B2")
        roc_df_out = roc_df.copy()
        roc_df_out["mode"] = "native"
        roc_df_out["auroc"] = roc_auc_val
        roc_data_rows.append(roc_df_out)

ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("ROC Curve")
ax.legend(loc="lower right")
fig.tight_layout()
fig.savefig(snakemake.output.roc_pdf, dpi=300)
plt.close(fig)

if roc_data_rows:
    pd.concat(roc_data_rows).to_csv(snakemake.output.roc_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["fpr", "tpr", "mode", "auroc"]).to_csv(snakemake.output.roc_data, sep="\t", index=False)

# --- Fig 2: PR curves ---
fig, ax = plt.subplots(figsize=(5, 5))
pr_data_rows = []

if best_score_col:
    _, pr_result = build_roc_pr_data(results_path, truth_path, best_score_col, is_pval, window)
    if pr_result:
        pr_df, pr_auc_val = pr_result
        ax.plot(pr_df["recall"], pr_df["precision"], label=f"Native (AUC={pr_auc_val:.3f})", color="#0072B2")
        pr_df_out = pr_df.copy()
        pr_df_out["mode"] = "native"
        pr_df_out["prauc"] = pr_auc_val
        pr_data_rows.append(pr_df_out)

ax.set_xlabel("Recall")
ax.set_ylabel("Precision")
ax.set_title("Precision-Recall Curve")
ax.legend(loc="lower left")
fig.tight_layout()
fig.savefig(snakemake.output.pr_pdf, dpi=300)
plt.close(fig)

if pr_data_rows:
    pd.concat(pr_data_rows).to_csv(snakemake.output.pr_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["recall", "precision", "mode", "prauc"]).to_csv(snakemake.output.pr_data, sep="\t", index=False)

# --- Fig 3: Score distribution by label ---
fig, ax = plt.subplots(figsize=(6, 4))
dist_data = pd.DataFrame()

if best_score_col:
    sep = "," if results_path.endswith(".csv") else "\t"
    try:
        preds = pd.read_csv(results_path, sep=sep)
        preds = normalize_columns(preds)
        truth = pd.read_csv(truth_path, sep="\t")
        truth = normalize_columns(truth)
        if "label" in truth.columns:
            truth_pos = truth[truth["label"] != "-"]
        else:
            truth_pos = truth

        truth_dict = {}
        for _, row in truth_pos.iterrows():
            tx = row["transcript"]
            if tx not in truth_dict:
                truth_dict[tx] = []
            truth_dict[tx].append(int(row["position"]))

        truth_txs = set(truth_pos["transcript"].unique())
        pred_on_truth = preds[preds["transcript"].isin(truth_txs)].copy()

        labels = []
        scores = []
        for _, row in pred_on_truth.iterrows():
            tx = row["transcript"]
            pos = int(row["position"])
            label = 0
            if tx in truth_dict:
                for tp in truth_dict[tx]:
                    if abs(pos - tp) <= window:
                        label = 1
                        break
            try:
                s = float(row[best_score_col])
            except (ValueError, TypeError):
                continue
            if np.isnan(s):
                continue
            labels.append(label)
            if is_pval:
                scores.append(-np.log10(max(s, 1e-300)))
            else:
                scores.append(s)

        if labels:
            dist_data = pd.DataFrame({"score": scores, "label": labels})
            pos_scores = dist_data[dist_data["label"] == 1]["score"]
            neg_scores = dist_data[dist_data["label"] == 0]["score"]
            if len(pos_scores) > 0:
                ax.hist(pos_scores, bins=30, alpha=0.6, label="Positive", color="#D55E00")
            if len(neg_scores) > 0:
                ax.hist(neg_scores, bins=30, alpha=0.6, label="Negative", color="#0072B2")
    except Exception:
        pass

score_label = f"-log10({best_score_col})" if is_pval else best_score_col
ax.set_xlabel(score_label if best_score_col else "Score")
ax.set_ylabel("Count")
ax.set_title("Score Distribution by Label")
ax.legend()
fig.tight_layout()
fig.savefig(snakemake.output.dist_pdf, dpi=300)
plt.close(fig)

if not dist_data.empty:
    dist_data.to_csv(snakemake.output.dist_data, sep="\t", index=False)
else:
    pd.DataFrame(columns=["score", "label"]).to_csv(snakemake.output.dist_data, sep="\t", index=False)

# --- Fig 4: Native vs Fair metric comparison ---
fig, ax = plt.subplots(figsize=(6, 4))
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

if comparison_data:
    comp_df = pd.DataFrame(comparison_data)
    comp_df.to_csv(snakemake.output.nvf_data, sep="\t", index=False)

    native_vals = comp_df[comp_df["mode"] == "native"].set_index("metric")["value"]
    fair_vals = comp_df[comp_df["mode"] == "fair"].set_index("metric")["value"]

    x = np.arange(len(metrics))
    width = 0.35
    if not native_vals.empty:
        ax.bar(x - width / 2, [native_vals.get(m, 0) for m in metrics], width, label="Native", color="#0072B2")
    if not fair_vals.empty:
        ax.bar(x + width / 2, [fair_vals.get(m, 0) for m in metrics], width, label="Fair", color="#D55E00")

    ax.set_xticks(x)
    ax.set_xticklabels([m.upper() for m in metrics], rotation=45, ha="right")
    ax.set_ylabel("Score")
    ax.set_title("Native vs Fair Mode Comparison")
    ax.legend()
    ax.set_ylim(0, 1.1)
else:
    pd.DataFrame(columns=["mode", "metric", "value"]).to_csv(snakemake.output.nvf_data, sep="\t", index=False)

fig.tight_layout()
fig.savefig(snakemake.output.nvf_pdf, dpi=300)
plt.close(fig)
