"""
Compute per-tool accuracy metrics against a ground-truth modification set.

Truth set TSV format (tab-separated, with header):
    transcript        - transcript identifier (must match tool output)
    position          - 0-based position on transcript
    modification_type - e.g. "m6A", "Psi", "m1Psi"
    label             - (optional) "-" for negative sites, any other value for positive
                       If absent, all sites are treated as positive.
                       Negative sites will be auto-inferred if not explicitly provided.

Tool result TSV format (standard NanoRNAMod output):
    transcript  - transcript identifier (or id, chrom, ref_id depending on tool)
    position    - 0-based position (or pos, start, start_loc depending on tool)
    (other tool-specific columns including optional scores)

Output columns:
    tool, modification_type, window, precision, recall, f1, tp, fp, fn, tn,
    specificity, mcc, auprc, auroc, called_sites,
    total_truth, total_predicted, total_negative

Configuration (from snakemake.params):
    window  - positional tolerance in nucleotides (0 = exact match)
               Can be a single integer or a list of integers for multi-window evaluation.
               Positions within ±window nt of a truth site count as a hit.
"""

import os
import ast
import pandas as pd
import numpy as np

result_files = snakemake.input.results
truth_set_path = snakemake.input.truth_set
output_files = snakemake.output  # Now expects 2 files
output_file = output_files[0]  # accuracy_summary.tsv (by modification_type)
output_overall = output_files[1]  # accuracy_summary_overall.tsv
window_param = snakemake.params.window

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
    # All sites are positive; negatives will be auto-inferred
    truth_pos = truth.copy()
    truth_neg = pd.DataFrame(columns=truth.columns)

# ── Helper: extract tool name from result file path ────────────────────────
def tool_from_path(path):
    """Extract tool name from path like .../modifications/{tool}/{comp}/{tool}_results.tsv"""
    parts = path.split(os.sep)
    for i, part in enumerate(parts):
        if part == "modifications" and i + 1 < len(parts):
            return parts[i + 1]
    return os.path.basename(os.path.dirname(os.path.dirname(path)))

# ── Helper: normalize column names for transcript and position ────────────
def normalize_columns(df):
    """Map various column names to standard 'transcript' and 'position'."""
    col_mapping = {}

    # Transcript column variants
    transcript_cols = ['transcript', 'id', 'ref_id', 'chrom']
    for col in transcript_cols:
        if col in df.columns:
            col_mapping[col] = 'transcript'
            break

    # Position column variants
    position_cols = ['position', 'pos', 'start', 'start_loc', 'transcript_pos', 'transcript_loc']
    for col in position_cols:
        if col in df.columns and col not in col_mapping.values():
            col_mapping[col] = 'position'
            break

    if col_mapping:
        df = df.rename(columns=col_mapping)

    return df

# ── Helper: detect score column in dataframe ───────────────────────────────
def detect_score_column(df, tool_name=None):
    """Detect the most likely score column for ranking predictions."""
    # Tool-specific score columns (in order of preference)
    tool_score_map = {
        'xpore': ['p_value', 'diff_mod', 'diff_mod_frac', 'mod_ratio'],
        'nanocompore': ['pvalue', 'logit_pvalue', 'logit', 'p_value', 'coverage'],
        'baleen': ['mod_score', 'score', 'kmer_score'],
        'differr': ['-log10 P value', '-log10_pvalue', 'score', 'pvalue'],
        'eligos2': ['pvalue', 'p_value', 'esb', 'oddsR'],
        'epinano': ['z_score_prediction', 'z_scores', 'delta_sum_err'],
        'drummer': ['p_value', 'pvalue', 'z_score'],
        'psipore': ['pvalue', 'p_value', 'score'],
        'tandemmod': ['probability', 'prob', 'score', 'mod_prob'],
        'directrm': ['probability', 'prob', 'mod_prob'],
        'm6atm': ['stoichiometry', 'probability', 'prob'],
        'rnano': ['probability', 'score', 'prob'],
        'nanopsu': ['pvalue', 'p_value', 'score'],
        'nanomud': ['probability', 'pvalue', 'score'],
        'penguin': ['probability', 'score', 'pvalue'],
    }

    # If tool name is provided, try tool-specific columns first
    if tool_name and tool_name in tool_score_map:
        for col in tool_score_map[tool_name]:
            # Try exact match
            if col in df.columns:
                return col
            # Try case-insensitive match
            for df_col in df.columns:
                if df_col.lower() == col.lower():
                    return df_col

    # Generic score columns (fallback)
    generic_score_cols = [
        'p_value', 'pvalue', 'pval', '-log10_pvalue', '-log10 P value',
        'score', 'probability', 'prob', 'mod_score', 'mod_prob',
        'logit_pvalue', 'logit', 'z_score', 'z_scores'
    ]

    for col in generic_score_cols:
        # Try exact match
        if col in df.columns:
            return col
        # Try case-insensitive match
        for df_col in df.columns:
            if df_col.lower().replace(' ', '_') == col.lower().replace(' ', '_'):
                return df_col

    return None

# ── Helper: compute AUPRC and AUROC ───────────────────────────────────────
def compute_ranking_metrics(pred_df, truth_pos_subset, truth_neg_subset, score_col, window=0):
    """
    Compute AUPRC and AUROC given predictions with scores.

    Uses binary labels based on match with ground truth sites within window.
    Score is converted to -log10(p-value) if it's a p-value column.
    """
    if score_col is None or pred_df.empty:
        return np.nan, np.nan

    # Check if score column exists and has valid values
    if score_col not in pred_df.columns:
        return np.nan, np.nan

    # Determine if this is a p-value column (needs -log10 transformation)
    score_col_lower = score_col.lower().replace(' ', '_')
    is_pvalue_col = any(x in score_col_lower for x in ['p_value', 'pvalue', 'pval'])

    # Create labels and scores arrays
    labels = []
    scores = []

    for _, pred_row in pred_df.iterrows():
        tx = pred_row['transcript']
        pos = pred_row['position']
        raw_score = pred_row[score_col]

        # Skip if score is NaN
        if pd.isna(raw_score):
            continue

        # Convert p-value to -log10(p-value) for ranking (higher = better)
        if is_pvalue_col:
            # Clamp p-value to avoid log(0)
            pval = max(float(raw_score), 1e-300)
            score = -np.log10(pval)
        else:
            score = float(raw_score)

        # Check if this prediction matches a positive truth site (within window)
        is_positive = False
        for _, truth_row in truth_pos_subset.iterrows():
            if (truth_row['transcript'] == tx and
                abs(truth_row['position'] - pos) <= window):
                is_positive = True
                break

        labels.append(1 if is_positive else 0)
        scores.append(score)

    # Need at least one positive and one negative sample
    if len(set(labels)) < 2 or len(scores) == 0:
        return np.nan, np.nan

    try:
        from sklearn.metrics import roc_auc_score, average_precision_score
        auroc = roc_auc_score(labels, scores)
        auprc = average_precision_score(labels, scores)
        return auprc, auroc
    except Exception as e:
        return np.nan, np.nan

# ── Load all tool results (preserving all columns for scores) ────────────────
tool_dfs = {}
tool_score_cols = {}

for f in result_files:
    tool = tool_from_path(f)
    try:
        df = pd.read_csv(f, sep='\t')
        df = normalize_columns(df)

        if 'transcript' in df.columns and 'position' in df.columns:
            # Detect score column for this tool (pass tool name for better detection)
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

# Get all modification types from positive sites
all_mod_types = set(truth_pos['modification_type'].unique())

for tool, pred_df in tool_dfs.items():
    for mod_type in all_mod_types:
        truth_subset = truth_pos[truth_pos['modification_type'] == mod_type].copy()

        if truth_subset.empty:
            continue

        # Get unique transcripts for this mod type
        all_tx = set(truth_subset['transcript'].unique())

        # Filter predictions to relevant transcripts
        pred_subset = pred_df[pred_df['transcript'].isin(all_tx)].copy()

        # Auto-infer negative sites if not provided
        # Negative sites = all positions that are NOT positive sites
        # (using the prediction set as the universe of possible positions)
        truth_neg_subset = truth_neg[truth_neg['modification_type'] == mod_type].copy()
        is_inferred_negative = False  # Track if negatives are explicit or inferred

        if truth_neg_subset.empty and not pred_subset.empty:
            # Create inferred negative sites from predictions that don't match truth
            neg_sites = []
            for _, pred_row in pred_subset.iterrows():
                tx = pred_row['transcript']
                pos = pred_row['position']
                is_match = False
                for _, truth_row in truth_subset.iterrows():
                    if (truth_row['transcript'] == tx and
                        abs(truth_row['position'] - pos) <= 0):
                        is_match = True
                        break
                if not is_match:
                    neg_sites.append({
                        'transcript': tx,
                        'position': pos,
                        'modification_type': mod_type,
                        'label': '-'
                    })
            if neg_sites:
                truth_neg_subset = pd.DataFrame(neg_sites)
                is_inferred_negative = True  # Mark as inferred

        score_col = tool_score_cols.get(tool)

        for window in windows:
            # Track which predicted sites matched truth sites
            matched_pred_indices = set()

            # For each positive truth site, check if there's a prediction within window
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

            called_sites = len(pred_subset)
            total_predicted = called_sites
            total_truth = len(truth_subset)

            # Calculate TN (true negatives)
            # TN = negative sites where there is NO prediction within window
            # Note: For inferred negatives, TN is not meaningful (negatives are derived from predictions)
            tn = 0
            if not truth_neg_subset.empty and not is_inferred_negative:
                # Only calculate TN for explicit negative sites
                for _, neg_row in truth_neg_subset.iterrows():
                    tx = neg_row['transcript']
                    pos = neg_row['position']

                    nearby_preds = pred_subset[
                        (pred_subset['transcript'] == tx) &
                        (pred_subset['position'].between(pos - window, pos + window))
                    ]

                    if nearby_preds.empty:
                        tn += 1

            total_negative = len(truth_neg_subset)

            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = (2 * precision * recall / (precision + recall)
                  if (precision + recall) > 0 else 0)

            # Specificity: TN / (TN + FP) or TN / total_negative
            # For inferred negatives, specificity is not meaningful
            if is_inferred_negative or total_negative == 0:
                specificity = np.nan
            else:
                specificity = tn / total_negative

            # Matthews Correlation Coefficient
            # For inferred negatives, MCC with TN=0 is misleading; set to NaN
            if is_inferred_negative:
                mcc = np.nan
            else:
                mcc = compute_mcc(tp, fp, fn, tn)

            # AUPRC and AUROC (only if score column exists)
            auprc, auroc = compute_ranking_metrics(
                pred_subset, truth_subset, truth_neg_subset, score_col, window
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
            })

# ── Write output: by modification_type ───────────────────────────────────────
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
        "total_truth", "total_predicted", "total_negative"
    ])
out_df.to_csv(output_file, sep='\t', index=False)

# ── Write output: overall (aggregated across all modification_types) ──────────
# Group by tool and window, summing TP/FP/FN/TN across modification types
overall_records = []
for tool in set(r["tool"] for r in records):
    for window in windows:
        tool_window_records = [r for r in records if r["tool"] == tool and r["window"] == window]
        if not tool_window_records:
            continue

        # Sum confusion matrix values
        tp_sum = sum(r["tp"] for r in tool_window_records)
        fp_sum = sum(r["fp"] for r in tool_window_records)
        fn_sum = sum(r["fn"] for r in tool_window_records)
        tn_sum = sum(r["tn"] for r in tool_window_records if not pd.isna(r["tn"]))
        tn_sum = tn_sum if tn_sum > 0 else 0

        total_truth_sum = sum(r["total_truth"] for r in tool_window_records)
        total_predicted_sum = sum(r["total_predicted"] for r in tool_window_records)
        total_negative_sum = sum(r["total_negative"] for r in tool_window_records)
        called_sites_sum = sum(r["called_sites"] for r in tool_window_records)

        # Check if any record has inferred negatives
        any_inferred = any(pd.isna(r.get("specificity")) for r in tool_window_records)

        # Calculate metrics
        precision = tp_sum / (tp_sum + fp_sum) if (tp_sum + fp_sum) > 0 else 0
        recall = tp_sum / (tp_sum + fn_sum) if (tp_sum + fn_sum) > 0 else 0
        f1 = (2 * precision * recall / (precision + recall)
              if (precision + recall) > 0 else 0)

        # Specificity and MCC
        if any_inferred or total_negative_sum == 0:
            specificity = np.nan
            mcc = np.nan
        else:
            specificity = tn_sum / total_negative_sum if total_negative_sum > 0 else np.nan
            mcc = compute_mcc(tp_sum, fp_sum, fn_sum, tn_sum)

        # AUPRC/AUROC: weighted average across modification types
        auprc_vals = [r["auprc"] for r in tool_window_records if not pd.isna(r["auprc"])]
        auroc_vals = [r["auroc"] for r in tool_window_records if not pd.isna(r["auroc"])]
        auprc_avg = np.mean(auprc_vals) if auprc_vals else np.nan
        auroc_avg = np.mean(auroc_vals) if auroc_vals else np.nan

        overall_records.append({
            "tool": tool,
            "window": window,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "tp": tp_sum, "fp": fp_sum, "fn": fn_sum, "tn": tn_sum,
            "specificity": specificity,
            "mcc": mcc,
            "auprc": auprc_avg,
            "auroc": auroc_avg,
            "called_sites": called_sites_sum,
            "total_truth": total_truth_sum,
            "total_predicted": total_predicted_sum,
            "total_negative": total_negative_sum,
        })

# Write overall output
if overall_records:
    overall_df = pd.DataFrame(overall_records).sort_values(
        ["window", "f1"],
        ascending=[True, False]
    )
else:
    overall_df = pd.DataFrame(columns=[
        "tool", "window", "precision", "recall",
        "f1", "tp", "fp", "fn", "tn", "specificity", "mcc",
        "auprc", "auroc", "called_sites",
        "total_truth", "total_predicted", "total_negative"
    ])
overall_df.to_csv(output_overall, sep='\t', index=False)

