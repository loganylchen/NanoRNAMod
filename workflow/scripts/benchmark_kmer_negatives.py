"""
K-mer Negative Benchmarking for RNA Modification Detection.

When ground truth does not contain explicit negative sites, we use k-mer matched
sites as negative controls. The rationale:
- RNA modifications often occur in specific sequence contexts (e.g., DRACH for m6A)
- Sites with the same k-mer context but NOT in the ground truth can serve as
  "biologically plausible" negative controls
- This enables meaningful MCC and specificity calculations

Input:
- Tool result TSVs (must contain 'kmer' or 'ref_kmer' column)
- Truth set TSV (positive sites only)

Output:
- accuracy_summary_kmer_negatives.tsv: Metrics using k-mer negatives
- kmer_negative_sites.tsv: The k-mer negative sites used for benchmarking
"""

import os
import ast
import pandas as pd
import numpy as np
from collections import defaultdict

# Snakemake inputs
result_files = snakemake.input.results
truth_set_path = snakemake.input.truth_set
output_metrics = snakemake.output.metrics
output_negatives = snakemake.output.negatives
window_param = snakemake.params.window
kmer_col_param = snakemake.params.get("kmer_column", "auto")

# Support both single int and list of windows
if isinstance(window_param, int):
    windows = [window_param]
elif isinstance(window_param, list):
    windows = window_param
else:
    try:
        parsed = ast.literal_eval(str(window_param))
        if isinstance(parsed, list):
            windows = parsed
        elif isinstance(parsed, int):
            windows = [parsed]
        else:
            windows = [int(window_param)]
    except (ValueError, SyntaxError):
        windows = [int(window_param)]


# Import shared utilities
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_utils import tool_from_path, normalize_columns


# ── Helper Functions ─────────────────────────────────────────────────────────

def detect_kmer_column(df, tool_name=None):
    """Detect the k-mer column in the dataframe."""
    kmer_cols = ['kmer', 'ref_kmer', 'k_mer', 'context', 'sequence_context', 'motif']

    for col in kmer_cols:
        if col in df.columns:
            return col

    df_cols_lower = {c.lower(): c for c in df.columns}
    for col in kmer_cols:
        if col.lower() in df_cols_lower:
            return df_cols_lower[col.lower()]

    return None


def detect_score_column(df, tool_name=None):
    """Detect the most likely score column for ranking predictions."""
    tool_score_patterns = {
        'xpore': ['pval_CASE_vs_CONTROL', 'z_score_CASE_vs_CONTROL', 'diff_mod_rate_CASE_vs_CONTROL',
                  'p_value', 'adjusted_p_value', 'statistic', 'pval', 'padj'],
        'nanocompore': ['GMM_logit_pvalue', 'KS_dwell_pvalue', 'KS_intensity_pvalue', 'Logit_LOR',
                        'p_value_glm', 'p_value_ks', 'p_value'],
        'baleen': ['p_value', 'padj', 'pvalue', 'effect_size', 'score'],
        'pybaleen': ['pvalue', 'padj', 'effect_size', 'p_value', 'adjusted_p_value', 'stoichiometry'],
        'differr': ['-log10 P value', '-log10 FDR', 'odds ratio', 'G statistic', 'score'],
        'drummer': ['OR_padj', 'G_padj', 'G_test', 'log2_(OR)', 'p_value', 'pvalue', 'padj'],
    }

    def match_column(df_columns, pattern):
        pattern_lower = pattern.lower().replace(' ', '_')
        if pattern in df_columns:
            return pattern
        for col in df_columns:
            if col.lower().replace(' ', '_') == pattern_lower:
                return col
        for col in df_columns:
            col_lower = col.lower().replace(' ', '_')
            if pattern_lower in col_lower:
                return col
        return None

    df_columns = list(df.columns)

    if tool_name and tool_name in tool_score_patterns:
        for pattern in tool_score_patterns[tool_name]:
            matched = match_column(df_columns, pattern)
            if matched:
                return matched

    generic_patterns = [
        'p_value', 'pvalue', 'pval', '-log10_pvalue', '-log10 P value',
        'score', 'probability', 'prob', 'mod_score', 'mod_prob',
        'logit_pvalue', 'logit', 'z_score', 'z_scores'
    ]

    for pattern in generic_patterns:
        matched = match_column(df_columns, pattern)
        if matched:
            return matched

    return None


def compute_ranking_metrics(pred_df, truth_pos_subset, truth_neg_subset, score_col, window=0, tool_name="unknown"):
    """Compute AUPRC and AUROC given predictions with scores."""
    if score_col is None or pred_df.empty:
        return np.nan, np.nan

    if score_col not in pred_df.columns:
        return np.nan, np.nan

    score_col_lower = score_col.lower().replace(' ', '_')
    is_log_transformed = '-log10' in score_col_lower or 'log10_pvalue' in score_col_lower
    is_pvalue_col = not is_log_transformed and any(x in score_col_lower for x in ['p_value', 'pvalue', 'pval', 'padj', 'q_value', 'qval', 'fdr', 'adj'])

    labels = []
    scores = []

    for _, pred_row in pred_df.iterrows():
        tx = pred_row['transcript']
        pos = pred_row['position']
        raw_score = pred_row[score_col]

        if pd.isna(raw_score):
            continue
        try:
            raw_score = float(raw_score)
        except (ValueError, TypeError):
            continue

        if is_pvalue_col:
            pval = max(raw_score, 1e-300)
            score = -np.log10(pval)
        else:
            score = raw_score

        is_positive = False
        for _, truth_row in truth_pos_subset.iterrows():
            if (truth_row['transcript'] == tx and
                abs(truth_row['position'] - pos) <= window):
                is_positive = True
                break

        is_negative = False
        if not is_positive and not truth_neg_subset.empty:
            for _, neg_row in truth_neg_subset.iterrows():
                if (neg_row['transcript'] == tx and
                    abs(neg_row['position'] - pos) <= window):
                    is_negative = True
                    break

        if is_positive:
            labels.append(1)
            scores.append(score)
        elif is_negative:
            labels.append(0)
            scores.append(score)

    if len(set(labels)) < 2 or len(scores) == 0:
        return np.nan, np.nan

    try:
        from sklearn.metrics import roc_auc_score, average_precision_score
        auroc = roc_auc_score(labels, scores)
        auprc = average_precision_score(labels, scores)
        return auprc, auroc
    except Exception:
        return np.nan, np.nan


def compute_mcc(tp, fp, fn, tn):
    """Compute Matthews Correlation Coefficient."""
    numerator = (tp * tn) - (fp * fn)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator == 0:
        return 0.0
    return numerator / denominator


# ── Load Truth Set ─────────────────────────────────────────────────────────
truth = pd.read_csv(truth_set_path, sep='\t')
required_cols = {"transcript", "position", "modification_type"}
assert required_cols.issubset(truth.columns), (
    f"Truth set must have columns: {required_cols}. Found: {set(truth.columns)}"
)

has_label = 'label' in truth.columns

if has_label:
    truth_neg_explicit = truth[truth['label'] == '-'].copy()
    if not truth_neg_explicit.empty:
        print(f"WARNING: Truth set has {len(truth_neg_explicit)} explicit negative sites.")
        print("K-mer negative benchmarking is most useful when ground truth lacks explicit negatives.")
    truth_pos = truth[truth['label'] != '-'].copy()
else:
    truth_pos = truth.copy()

print(f"Loaded {len(truth_pos)} positive sites from truth set")


# ── Load All Tool Results ─────────────────────────────────────────────────────
tool_dfs = {}
tool_kmer_cols = {}
tool_score_cols = {}

for f in result_files:
    tool = tool_from_path(f)
    try:
        df = pd.read_csv(f, sep='\t')
        df = normalize_columns(df)

        if 'transcript' in df.columns and 'position' in df.columns:
            kmer_col = detect_kmer_column(df, tool)
            if kmer_col and tool not in tool_kmer_cols:
                tool_kmer_cols[tool] = kmer_col
                print(f"Detected k-mer column for {tool}: {kmer_col}")

            score_col = detect_score_column(df, tool)
            if score_col and tool not in tool_score_cols:
                tool_score_cols[tool] = score_col
                print(f"Detected score column for {tool}: {score_col}")

            if tool not in tool_dfs:
                tool_dfs[tool] = []
            tool_dfs[tool].append(df)
    except Exception as e:
        print(f"Warning: could not read {f}: {e}")

for tool in tool_dfs:
    tool_dfs[tool] = pd.concat(tool_dfs[tool], ignore_index=True).drop_duplicates(
        subset=['transcript', 'position']
    )
    print(f"Loaded {len(tool_dfs[tool])} predictions for {tool}")


# ── Generate K-mer Negative Sites ─────────────────────────────────────────────
def generate_kmer_negatives(truth_pos, tool_dfs, tool_kmer_cols):
    """Generate k-mer matched negative sites from tool predictions."""
    kmer_negatives = defaultdict(list)
    positive_kmers_by_mod = defaultdict(set)

    for mod_type in truth_pos['modification_type'].unique():
        truth_subset = truth_pos[truth_pos['modification_type'] == mod_type]

        for _, row in truth_subset.iterrows():
            tx = row['transcript']
            pos = row['position']

            for tool, df in tool_dfs.items():
                if tool not in tool_kmer_cols:
                    continue
                kmer_col = tool_kmer_cols[tool]

                matches = df[
                    (df['transcript'] == tx) &
                    (df['position'] == pos)
                ]

                if not matches.empty:
                    kmer = matches.iloc[0][kmer_col]
                    if pd.notna(kmer):
                        positive_kmers_by_mod[(tx, mod_type)].add(str(kmer).upper())

    positive_sites = set()
    for _, row in truth_pos.iterrows():
        positive_sites.add((row['transcript'], row['position']))

    for mod_type in truth_pos['modification_type'].unique():
        all_kmers_for_mod = set()
        for (tx, mt), kmers in positive_kmers_by_mod.items():
            if mt == mod_type:
                all_kmers_for_mod.update(kmers)

        if not all_kmers_for_mod:
            continue

        for tool, df in tool_dfs.items():
            if tool not in tool_kmer_cols:
                continue
            kmer_col = tool_kmer_cols[tool]

            relevant_txs = set(truth_pos[truth_pos['modification_type'] == mod_type]['transcript'].unique())
            df_filtered = df[df['transcript'].isin(relevant_txs)]

            for _, row in df_filtered.iterrows():
                tx = row['transcript']
                pos = row['position']
                kmer = row[kmer_col]

                if pd.isna(kmer):
                    continue

                kmer_upper = str(kmer).upper()
                if kmer_upper in all_kmers_for_mod:
                    if (tx, pos) not in positive_sites:
                        kmer_negatives[mod_type].append({
                            'transcript': tx,
                            'position': pos,
                            'modification_type': mod_type,
                            'kmer': kmer,
                            'source_tool': tool,
                            'label': '-'
                        })

    all_negatives = []
    for mod_type, sites in kmer_negatives.items():
        seen = set()
        for site in sites:
            key = (site['transcript'], site['position'])
            if key not in seen:
                seen.add(key)
                all_negatives.append(site)

    return pd.DataFrame(all_negatives) if all_negatives else pd.DataFrame()


print("\nGenerating k-mer matched negative sites...")
kmer_negatives_df = generate_kmer_negatives(truth_pos, tool_dfs, tool_kmer_cols)

if kmer_negatives_df.empty:
    print("WARNING: Could not generate k-mer negative sites.")
    print("This usually means tool outputs don't contain k-mer information.")
    pd.DataFrame().to_csv(output_metrics, sep='\t', index=False)
    pd.DataFrame().to_csv(output_negatives, sep='\t', index=False)
    exit(0)

print(f"Generated {len(kmer_negatives_df)} k-mer matched negative sites")
for mod_type in kmer_negatives_df['modification_type'].unique():
    n = len(kmer_negatives_df[kmer_negatives_df['modification_type'] == mod_type])
    print(f"  {mod_type}: {n} negative sites")


# ── Compute Metrics with K-mer Negatives ───────────────────────────────────────
records = []
all_mod_types = set(truth_pos['modification_type'].unique())

for tool, pred_df in tool_dfs.items():
    for mod_type in all_mod_types:
        truth_subset = truth_pos[truth_pos['modification_type'] == mod_type].copy()

        if truth_subset.empty:
            continue

        kmer_neg_subset = kmer_negatives_df[
            kmer_negatives_df['modification_type'] == mod_type
        ].copy() if not kmer_negatives_df.empty else pd.DataFrame()

        all_tx = set(truth_subset['transcript'].unique())
        if not kmer_neg_subset.empty:
            all_tx.update(kmer_neg_subset['transcript'].unique())

        pred_subset = pred_df[pred_df['transcript'].isin(all_tx)].copy()

        score_col = tool_score_cols.get(tool)

        for window in windows:
            matched_pred_indices = set()

            tp = 0
            for _, truth_row in truth_subset.iterrows():
                tx = truth_row['transcript']
                pos = truth_row['position']

                matches = pred_subset[
                    (pred_subset['transcript'] == tx) &
                    (pred_subset['position'].between(pos - window, pos + window))
                ]

                if not matches.empty:
                    tp += 1
                    matched_pred_indices.update(matches.index.tolist())

            fn = len(truth_subset) - tp
            fp = len(pred_subset) - len(matched_pred_indices)

            tn = 0
            if not kmer_neg_subset.empty:
                for _, neg_row in kmer_neg_subset.iterrows():
                    tx = neg_row['transcript']
                    pos = neg_row['position']

                    nearby_preds = pred_subset[
                        (pred_subset['transcript'] == tx) &
                        (pred_subset['position'].between(pos - window, pos + window))
                    ]

                    if nearby_preds.empty:
                        tn += 1

            called_sites = len(pred_subset)
            total_predicted = called_sites
            total_truth = len(truth_subset)
            total_negative = len(kmer_neg_subset)

            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = (2 * precision * recall / (precision + recall)
                  if (precision + recall) > 0 else 0)

            specificity = tn / total_negative if total_negative > 0 else np.nan

            if total_negative > 0:
                mcc = compute_mcc(tp, fp, fn, tn)
            else:
                mcc = np.nan

            auprc, auroc = compute_ranking_metrics(
                pred_subset, truth_subset, kmer_neg_subset, score_col, window, tool_name=tool
            )

            records.append({
                "tool": tool,
                "modification_type": mod_type,
                "window": window,
                "precision": precision,
                "recall": recall,
                "f1": f1,
                "tp": tp, "fp": fp, "fn": fn, "tn": tn,
                "specificity": specificity,
                "mcc": mcc,
                "auprc": auprc,
                "auroc": auroc,
                "called_sites": called_sites,
                "total_truth": total_truth,
                "total_predicted": total_predicted,
                "total_negative": total_negative,
                "negative_type": "kmer_matched",
            })


# ── Write Output ─────────────────────────────────────────────────────────────
if records:
    out_df = pd.DataFrame(records).sort_values(
        ["modification_type", "window", "f1"],
        ascending=[True, True, False]
    )
else:
    out_df = pd.DataFrame(columns=[
        "tool", "modification_type", "window", "precision", "recall",
        "f1", "tp", "fp", "fn", "tn", "specificity", "mcc",
        "auprc", "auroc", "called_sites",
        "total_truth", "total_predicted", "total_negative", "negative_type"
    ])

out_df.to_csv(output_metrics, sep='\t', index=False)
print(f"\nSaved metrics to {output_metrics}")

kmer_negatives_df.to_csv(output_negatives, sep='\t', index=False)
print(f"Saved {len(kmer_negatives_df)} k-mer negative sites to {output_negatives}")
