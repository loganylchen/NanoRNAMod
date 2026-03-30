"""
Generate per-sample tool figures (native only, no fair mode comparison).

Produces a score distribution figure with co-located source data.

Inputs (via Snakemake):
    native_scores: score_comparison.tsv from native mode
    results: raw tool results TSV
    truth_set: truth set TSV

Outputs:
    dist_pdf: score_distribution.pdf
    dist_data: score_distribution_data.tsv
"""

import os
import pandas as pd
import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("matplotlib is required for per-sample figures")


def normalize_columns(df):
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


# ===========================================================================
# Main
# ===========================================================================

native_scores_path = snakemake.input.native_scores
results_path = snakemake.input.results
truth_path = snakemake.input.truth_set

window_param = snakemake.params.window
if isinstance(window_param, list):
    window = max(int(w) for w in window_param)
else:
    window = int(window_param)

out_dir = os.path.dirname(snakemake.output.dist_pdf)
os.makedirs(out_dir, exist_ok=True)

# Load native scores to find best score column
try:
    native_df = pd.read_csv(native_scores_path, sep="\t")
except Exception:
    native_df = pd.DataFrame()

best_score_col = None
is_pval = False
if not native_df.empty:
    best_score_col = native_df.iloc[0].get("score_column", None)
    is_pval = bool(native_df.iloc[0].get("is_pvalue", False))

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
        preds["position"] = pd.to_numeric(preds["position"], errors="coerce")
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
