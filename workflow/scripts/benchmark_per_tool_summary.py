"""
Generate per-tool summary figures across all comparisons.

For each tool, produces 4 figures with co-located source data:
  - accuracy_by_comparison.pdf — grouped bar of AUROC/F1 per comparison
  - native_vs_fair_summary.pdf — paired dot plot: native vs fair per comparison
  - score_columns_comparison.pdf — all score columns ranked
  - threshold_curve.pdf — F1 vs threshold for best score column

Inputs (via Snakemake):
    native_scores: list of score_comparison.tsv (one per comparison, native mode)
    fair_scores: list of score_comparison.tsv (one per comparison, fair mode)

Outputs:
    8 files in per_tool/{tool}/figures/ (4 PDFs + 4 data TSVs)
"""

import os
import pandas as pd
import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("matplotlib is required for per-tool summary figures")


def load_and_tag(paths, mode):
    """Load score_comparison.tsv files and tag with comparison + mode."""
    records = []
    for path in paths:
        parts = path.replace("\\", "/").split("/")
        # Extract comparison from: per_tool/{tool}/{comparison}/{mode}/score_comparison.tsv
        comp = None
        for i, part in enumerate(parts):
            if part == "per_tool" and i + 2 < len(parts):
                comp = parts[i + 2]
                break
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            continue
        if df.empty:
            continue
        df["comparison"] = comp
        df["mode"] = mode
        records.append(df)
    if records:
        return pd.concat(records, ignore_index=True)
    return pd.DataFrame()


def extract_comparison(path):
    """Extract comparison from per_tool/{tool}/{comp}/{mode}/score_comparison.tsv"""
    parts = path.replace("\\", "/").split("/")
    for i, part in enumerate(parts):
        if part == "per_tool" and i + 2 < len(parts):
            return parts[i + 2]
    return os.path.basename(os.path.dirname(os.path.dirname(path)))


# ===========================================================================
# Main
# ===========================================================================

native_paths = list(snakemake.input.native_scores)
fair_paths = list(snakemake.input.fair_scores)

out_dir = os.path.dirname(snakemake.output.acc_pdf)
os.makedirs(out_dir, exist_ok=True)

native_df = load_and_tag(native_paths, "native")
fair_df = load_and_tag(fair_paths, "fair")

all_df = pd.concat([native_df, fair_df], ignore_index=True) if not native_df.empty or not fair_df.empty else pd.DataFrame()

# --- Fig 1: Accuracy by Comparison (grouped bar) ---
fig, ax = plt.subplots(figsize=(8, 5))
acc_data = pd.DataFrame()

if not native_df.empty:
    # Use best score column per comparison (first row = highest AUROC)
    best_per_comp = native_df.groupby("comparison").first().reset_index()
    acc_data = best_per_comp[["comparison", "auroc", "f1"]].copy()
    acc_data["mode"] = "native"

    if not fair_df.empty:
        fair_best = fair_df.groupby("comparison").first().reset_index()
        fair_acc = fair_best[["comparison", "auroc", "f1"]].copy()
        fair_acc["mode"] = "fair"
        acc_data = pd.concat([acc_data, fair_acc], ignore_index=True)

if not acc_data.empty:
    comparisons = acc_data["comparison"].unique()
    x = np.arange(len(comparisons))
    width = 0.35

    native_auroc = acc_data[acc_data["mode"] == "native"].set_index("comparison")
    fair_auroc = acc_data[acc_data["mode"] == "fair"].set_index("comparison") if "fair" in acc_data["mode"].values else pd.DataFrame()

    ax.bar(x - width / 2, [native_auroc.loc[c, "auroc"] if c in native_auroc.index else 0 for c in comparisons],
           width, label="Native AUROC", color="#0072B2")
    if not fair_auroc.empty:
        ax.bar(x + width / 2, [fair_auroc.loc[c, "auroc"] if c in fair_auroc.index else 0 for c in comparisons],
               width, label="Fair AUROC", color="#D55E00")

    ax.set_xticks(x)
    ax.set_xticklabels(comparisons, rotation=45, ha="right")
    ax.set_ylabel("AUROC")
    ax.set_title("Accuracy by Comparison")
    ax.legend()
    ax.set_ylim(0, 1.1)

acc_data.to_csv(snakemake.output.acc_data, sep="\t", index=False)
fig.tight_layout()
fig.savefig(snakemake.output.acc_pdf, dpi=300)
plt.close(fig)

# --- Fig 2: Native vs Fair Summary (paired dot plot) ---
fig, ax = plt.subplots(figsize=(6, 5))
nvf_data = pd.DataFrame()

if not acc_data.empty and "fair" in acc_data["mode"].values:
    comparisons = sorted(acc_data["comparison"].unique())
    for i, comp in enumerate(comparisons):
        native_row = acc_data[(acc_data["comparison"] == comp) & (acc_data["mode"] == "native")]
        fair_row = acc_data[(acc_data["comparison"] == comp) & (acc_data["mode"] == "fair")]
        if not native_row.empty and not fair_row.empty:
            n_val = float(native_row.iloc[0]["auroc"])
            f_val = float(fair_row.iloc[0]["auroc"])
            ax.plot([i, i], [n_val, f_val], "k-", alpha=0.3)
            ax.scatter(i, n_val, color="#0072B2", s=60, zorder=5, label="Native" if i == 0 else "")
            ax.scatter(i, f_val, color="#D55E00", s=60, zorder=5, label="Fair" if i == 0 else "")

    ax.set_xticks(range(len(comparisons)))
    ax.set_xticklabels(comparisons, rotation=45, ha="right")
    ax.set_ylabel("AUROC")
    ax.set_title("Native vs Fair per Comparison")
    ax.legend()
    ax.set_ylim(0, 1.1)
    nvf_data = acc_data

nvf_data.to_csv(snakemake.output.nvf_data, sep="\t", index=False)
fig.tight_layout()
fig.savefig(snakemake.output.nvf_pdf, dpi=300)
plt.close(fig)

# --- Fig 3: Score Columns Comparison ---
fig, ax = plt.subplots(figsize=(8, 5))
sc_data = pd.DataFrame()

if not native_df.empty:
    # Aggregate mean AUROC per score column across comparisons
    sc_agg = (
        native_df.groupby("score_column")
        .agg(mean_auroc=("auroc", "mean"), std_auroc=("auroc", "std"), n=("comparison", "nunique"))
        .reset_index()
        .sort_values("mean_auroc", ascending=False)
    )
    sc_data = sc_agg.copy()

    ax.barh(range(len(sc_agg)), sc_agg["mean_auroc"], xerr=sc_agg["std_auroc"].fillna(0),
            color="#4477AA", alpha=0.8)
    ax.set_yticks(range(len(sc_agg)))
    ax.set_yticklabels(sc_agg["score_column"], fontsize=8)
    ax.set_xlabel("Mean AUROC")
    ax.set_title("Score Column Ranking")
    ax.set_xlim(0, 1.1)
    ax.invert_yaxis()

sc_data.to_csv(snakemake.output.sc_data, sep="\t", index=False)
fig.tight_layout()
fig.savefig(snakemake.output.sc_pdf, dpi=300)
plt.close(fig)

# --- Fig 4: Threshold Curve (F1 vs threshold for best score) ---
fig, ax = plt.subplots(figsize=(6, 4))
thresh_data = pd.DataFrame()

if not native_df.empty:
    best_col = native_df.iloc[0].get("score_column", None)
    if best_col:
        # Gather threshold/F1 across comparisons
        rows = []
        for _, row in native_df[native_df["score_column"] == best_col].iterrows():
            rows.append({
                "comparison": row.get("comparison", ""),
                "best_threshold": row.get("best_threshold", np.nan),
                "f1": row.get("f1", np.nan),
            })
        if rows:
            thresh_data = pd.DataFrame(rows)
            ax.scatter(thresh_data["best_threshold"], thresh_data["f1"],
                       color="#0072B2", s=60, zorder=5)
            for _, r in thresh_data.iterrows():
                ax.annotate(r["comparison"], (r["best_threshold"], r["f1"]),
                            fontsize=7, ha="left", va="bottom")

ax.set_xlabel("Threshold")
ax.set_ylabel("F1 Score")
ax.set_title("F1 vs Threshold (Best Score Column)")
thresh_data.to_csv(snakemake.output.thresh_data, sep="\t", index=False)
fig.tight_layout()
fig.savefig(snakemake.output.thresh_pdf, dpi=300)
plt.close(fig)
