"""
Compute per-tool accuracy metrics against a ground-truth modification set.

Truth set TSV format (tab-separated, with header):
    transcript        - transcript identifier (must match tool output)
    position          - 0-based position on transcript
    modification_type - e.g. "m6A", "Psi", "m1Psi"
    label             - (optional) "-" for negative sites, any other value for positive
                       If absent, all sites are treated as positive.

Tool result TSV format (standard NanoRNAMod output):
    transcript  - transcript identifier
    position    - 0-based position
    (other tool-specific columns)

Output columns:
    tool, modification_type, window, precision, recall, f1, tp, fp, fn, tn, specificity, mcc, total_truth, total_predicted, total_negative

Configuration (from snakemake.params):
    window  - positional tolerance in nucleotides (0 = exact match)
               Can be a single integer or a list of integers for multi-window evaluation.
               Positions within ±window nt of a truth site count as a hit.
"""

import os
import ast
import pandas as pd
import numpy as np

result_files = snakemake.input.results    # list of *_results.tsv paths
truth_set_path = snakemake.input.truth_set
output_file = snakemake.output[0]
window_param = snakemake.params.window

# Support both single int and list of windows
if isinstance(window_param, int):
    windows = [window_param]
elif isinstance(window_param, list):
    windows = window_param
else:
    # Try to parse as string representation (safe parsing with ast.literal_eval)
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

# ── Load truth set ─────────────────────────────────────────────────────────
truth = pd.read_csv(truth_set_path, sep='\t')
required_cols = {"transcript", "position", "modification_type"}
assert required_cols.issubset(truth.columns), (
    f"Truth set must have columns: {required_cols}. Found: {set(truth.columns)}"
)

# Check for label column (optional)
has_label = 'label' in truth.columns

if has_label:
    # Split truth into positive and negative sites
    truth_pos = truth[truth['label'] != '-'].copy()
    truth_neg = truth[truth['label'] == '-'].copy()
else:
    # All sites are positive
    truth_pos = truth.copy()
    truth_neg = pd.DataFrame(columns=truth.columns)

# ── Helper: extract tool name from result file path ────────────────────────
def tool_from_path(path):
    """
    Path pattern: .../modifications/{tool}/{comp_or_sample}/{tool}_results.tsv
    """
    parts = path.split(os.sep)
    for i, part in enumerate(parts):
        if part == "modifications" and i + 1 < len(parts):
            return parts[i + 1]
    return os.path.basename(os.path.dirname(os.path.dirname(path)))

# ── Load all tool results ──────────────────────────────────────────────────
tool_dfs = {}
for f in result_files:
    tool = tool_from_path(f)
    try:
        df = pd.read_csv(f, sep='\t')
        if 'transcript' in df.columns and 'position' in df.columns:
            if tool not in tool_dfs:
                tool_dfs[tool] = []
            tool_dfs[tool].append(df[['transcript', 'position']])
    except Exception as e:
        print(f"Warning: could not read {f}: {e}")

for tool in tool_dfs:
    tool_dfs[tool] = pd.concat(tool_dfs[tool], ignore_index=True).drop_duplicates()

# ─────────────────────────────────────────────────────────────────────────────
# Compute accuracy metrics for each tool, modification type, and window
# ─────────────────────────────────────────────────────────────────────────────

def compute_mcc(tp, fp, fn, tn):
    """Compute Matthews Correlation Coefficient."""
    numerator = (tp * tn) - (fp * fn)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator == 0:
        return 0.0
    return numerator / denominator

records = []

# Get all modification types from both positive and negative sites
all_mod_types = set(truth_pos['modification_type'].unique()) | set(truth_neg['modification_type'].unique())

for tool, pred_df in tool_dfs.items():
    for mod_type in all_mod_types:
        truth_subset = truth_pos[truth_pos['modification_type'] == mod_type].copy()
        truth_neg_subset = truth_neg[truth_neg['modification_type'] == mod_type].copy()

        if truth_subset.empty and truth_neg_subset.empty:
            continue

        # Get unique transcripts for this mod type
        tx_with_pos = set(truth_subset['transcript'].unique()) if not truth_subset.empty else set()
        tx_with_neg = set(truth_neg_subset['transcript'].unique()) if not truth_neg_subset.empty else set()
        all_tx = tx_with_pos | tx_with_neg

        # Filter predictions to relevant transcripts
        pred_subset = pred_df[pred_df['transcript'].isin(all_tx)].copy()

        for window in windows:
            # Track which predicted sites matched truth sites
            matched_pred_indices = set()

            # For each positive truth site, check if there's a prediction within window
            tp = 0
            for _, truth_row in truth_subset.iterrows():
                tx = truth_row['transcript']
                pos = truth_row['position']

                # Find predictions within window on same transcript
                matches = pred_subset[
                    (pred_subset['transcript'] == tx) &
                    (pred_subset['position'].between(pos - window, pos + window))
                ]

                if not matches.empty:
                    tp += 1
                    matched_pred_indices.update(matches.index.tolist())

            fn = len(truth_subset) - tp
            fp = len(pred_subset) - len(matched_pred_indices)

            total_predicted = len(pred_subset)
            total_truth = len(truth_subset)

            # Calculate TN (true negatives)
            # TN = negative sites where there is NO prediction within window
            tn = 0
            if not truth_neg_subset.empty:
                for _, neg_row in truth_neg_subset.iterrows():
                    tx = neg_row['transcript']
                    pos = neg_row['position']

                    # Check if there are predictions within window of this negative site
                    nearby_preds = pred_subset[
                        (pred_subset['transcript'] == tx) &
                        (pred_subset['position'].between(pos - window, pos + window))
                    ]

                    # If no predictions nearby, this is a true negative
                    if nearby_preds.empty:
                        tn += 1

            total_negative = len(truth_neg_subset)

            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = (2 * precision * recall / (precision + recall)
                  if (precision + recall) > 0 else 0)

            # Specificity: TN / (TN + FP)
            # Note: FP here includes all false positives, not just those near negative sites
            if total_negative > 0:
                specificity = tn / total_negative
            else:
                specificity = np.nan

            # Matthews Correlation Coefficient
            mcc = compute_mcc(tp, fp, fn, tn)

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
                "total_truth": total_truth,
                "total_predicted": total_predicted,
                "total_negative": total_negative,
            })

# Write output
if records:
    out_df = pd.DataFrame(records).sort_values(["modification_type", "window", "f1"], ascending=[True, True, False])
else:
    out_df = pd.DataFrame(columns=["tool", "modification_type", "window", "precision", "recall",
                                   "f1", "tp", "fp", "fn", "tn", "specificity", "mcc",
                                   "total_truth", "total_predicted", "total_negative"])
out_df.to_csv(output_file, sep='\t', index=False)
