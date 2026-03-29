"""
Same-Base Negative Benchmarking for RNA Modification Detection.

When ground truth does not contain explicit negative sites, we use same-base
sites as negative controls. The rationale:
- RNA modifications occur at specific nucleotides (e.g., m6A at adenines)
- Sites with the same base but NOT in the ground truth serve as
  "sequence-matched" negative controls
- This is Strategy 2: less stringent than k-mer matching but more stringent
  than using all non-positive sites
- This enables meaningful MCC and specificity calculations

Input:
- Tool result TSVs (must contain 'ref_base', 'base', 'nucleotide', or similar column)
- Truth set TSV (positive sites only)

Output:
- accuracy_summary_same_base_negatives.tsv: Metrics using same-base negatives
- same_base_negative_sites.tsv: The same-base negative sites used for benchmarking
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
base_col_param = snakemake.params.get("base_column", "auto")

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




def tool_from_path(path):
    """Extract tool name from file path."""
    parts = path.split(os.sep)
    for i, part in enumerate(parts):
        if part == "modifications" and i + 1 < len(parts):
            return parts[i + 1]
    return os.path.basename(os.path.dirname(os.path.dirname(path)))


def normalize_columns(df):
    """Normalize various column names to standard 'transcript' and 'position'."""
    col_mapping = {}
    transcript_cols = ['transcript', 'id', 'ref_id', 'chrom', 'transcript_id']
    for col in transcript_cols:
        if col in df.columns:
            col_mapping[col] = 'transcript'
            break
    position_cols = ['position', 'pos', 'start', 'start_loc', 'transcript_pos', 'transcript_loc']
    for col in position_cols:
        if col in df.columns and col not in col_mapping.values():
            col_mapping[col] = 'position'
            break
    if col_mapping:
        df = df.rename(columns=col_mapping)
    return df


# ── Helper Functions ─────────────────────────────────────────────────────────

def detect_base_column(df, tool_name=None):
    """Detect the base/nucleotide column in the dataframe."""
    base_cols = ['ref_base', 'base', 'nucleotide', 'ref_nucleotide', 'nt', 'ref_nt',
                 'modified_base', 'mod_base']

    for col in base_cols:
        if col in df.columns:
            return col

    df_cols_lower = {c.lower(): c for c in df.columns}
    for col in base_cols:
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
        'pybaleen': ['mod_ratio', 'pvalue', 'padj', 'effect_size', 'mean_p_mod',
                      'p_value', 'adjusted_p_value', 'stoichiometry'],
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

        # Check if positive
        is_positive = False
        for _, truth_row in truth_pos_subset.iterrows():
            if (truth_row['transcript'] == tx and
                abs(truth_row['position'] - pos) <= window):
                is_positive = True
                break

        # Check if negative
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
    if tp + fp == 0 or tp + fn == 0 or tn + fp == 0 or tn + fn == 0:
        return 0.0
    numerator = (tp * tn) - (fp * fn)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator == 0:
        return 0.0
    return numerator / denominator


# ── Load Ground Truth ────────────────────────────────────────────────────────

truth_df = pd.read_csv(truth_set_path, sep='\t')
truth_df = normalize_columns(truth_df)

if 'modification_type' not in truth_df.columns:
    truth_df['modification_type'] = 'unknown'

# Split into positive and negative (if provided)
if 'site_type' in truth_df.columns:
    truth_pos = truth_df[truth_df['site_type'].isin(['positive', 'Positive', 'POS', 'pos', 1, '1'])].copy()
    truth_neg = truth_df[truth_df['site_type'].isin(['negative', 'Negative', 'NEG', 'neg', 0, '0'])].copy()
else:
    truth_pos = truth_df.copy()
    truth_neg = pd.DataFrame()

print(f"Loaded {len(truth_pos)} positive truth sites")
print(f"Loaded {len(truth_neg)} explicit negative truth sites")


# ── Load Tool Results ────────────────────────────────────────────────────────

tool_dfs = {}
tool_score_cols = {}
tool_base_cols = {}

for f in result_files:
    tool = tool_from_path(f)
    try:
        df = pd.read_csv(f, sep='\t')
        df = normalize_columns(df)

        if 'transcript' in df.columns and 'position' in df.columns:
            # Detect columns
            base_col = detect_base_column(df, tool)
            score_col = detect_score_column(df, tool)

            if base_col:
                df['_base_col'] = base_col

            if tool not in tool_dfs:
                tool_dfs[tool] = []
            tool_dfs[tool].append(df)

            if score_col and tool not in tool_score_cols:
                tool_score_cols[tool] = score_col

            if base_col and tool not in tool_base_cols:
                tool_base_cols[tool] = base_col

            print(f"Loaded {tool}: {len(df)} predictions, base_col={base_col}, score_col={score_col}")
    except Exception as e:
        print(f"Warning: could not read {f}: {e}")

for tool in tool_dfs:
    tool_dfs[tool] = pd.concat(tool_dfs[tool], ignore_index=True).drop_duplicates(
        subset=['transcript', 'position']
    )


# ── Generate Same-Base Negative Sites ────────────────────────────────────────
# Strategy: For each modification type, collect all sites with the same base
# as positive sites but NOT in the positive set

def generate_same_base_negatives(tool_dfs, truth_pos, tool_base_cols):
    """
    Generate negative sites by finding positions with the same base as positive sites
    but not in the positive set.

    For m6A: find all adenines (A) in tool predictions that are not positive sites
    """
    all_negatives = []
    all_mod_types = truth_pos['modification_type'].unique()

    for mod_type in all_mod_types:
        mod_truth = truth_pos[truth_pos['modification_type'] == mod_type]

        if mod_truth.empty:
            continue

        # Collect positive site positions
        positive_sites = set()
        for _, row in mod_truth.iterrows():
            positive_sites.add((row['transcript'], row['position']))

        # Infer the modified base from positive sites
        # First try to get from truth set
        modified_bases = set()
        if 'ref_base' in mod_truth.columns:
            modified_bases.update(mod_truth['ref_base'].dropna().unique())
        elif 'base' in mod_truth.columns:
            modified_bases.update(mod_truth['base'].dropna().unique())

        # If no base info in truth, infer from tool results
        if not modified_bases:
            for tool, df in tool_dfs.items():
                base_col = tool_base_cols.get(tool)
                if base_col and base_col in df.columns:
                    # Find bases at positive sites
                    for _, row in mod_truth.iterrows():
                        tx = row['transcript']
                        pos = row['position']
                        matches = df[(df['transcript'] == tx) & (df['position'] == pos)]
                        if not matches.empty:
                            bases = matches[base_col].dropna().unique()
                            modified_bases.update(bases)

        # If still no base info, use common modification bases
        if not modified_bases:
            # Default bases for common modifications
            mod_type_lower = str(mod_type).lower()
            if 'm6a' in mod_type_lower or 'm1a' in mod_type_lower:
                modified_bases = {'A', 'a'}
            elif 'm5c' in mod_type_lower or 'm3c' in mod_type_lower or 'hm5c' in mod_type_lower:
                modified_bases = {'C', 'c'}
            elif 'psu' in mod_type_lower or 'pseudouridine' in mod_type_lower:
                modified_bases = {'U', 'u', 'T', 't'}
            elif 'm7g' in mod_type_lower:
                modified_bases = {'G', 'g'}
            else:
                modified_bases = {'A', 'a', 'C', 'c', 'G', 'g', 'U', 'u', 'T', 't'}

        modified_bases = {b.upper() for b in modified_bases if pd.notna(b)}
        print(f"  {mod_type}: Modified bases = {modified_bases}")

        # Find candidate negative sites from tool predictions
        # These are sites with the same base but not in positive set
        candidate_sites = {}  # (tx, pos) -> base

        for tool, df in tool_dfs.items():
            base_col = tool_base_cols.get(tool)
            if not base_col or base_col not in df.columns:
                continue

            # Filter to same-base sites
            base_mask = df[base_col].apply(lambda x: str(x).upper() in modified_bases if pd.notna(x) else False)
            same_base_df = df[base_mask]

            for _, row in same_base_df.iterrows():
                tx = row['transcript']
                pos = row['position']
                site_key = (tx, pos)

                # Skip if in positive set
                if site_key in positive_sites:
                    continue

                candidate_sites[site_key] = row[base_col]

        # Convert to negative sites list
        for (tx, pos), base in candidate_sites.items():
            all_negatives.append({
                'transcript': tx,
                'position': pos,
                'modification_type': mod_type,
                'ref_base': base,
                'site_type': 'same_base_negative'
            })

        print(f"  {mod_type}: Generated {len(candidate_sites)} same-base negative sites")

    return pd.DataFrame(all_negatives)


print("\nGenerating same-base negative sites...")
if not truth_neg.empty:
    print("Using explicit negative sites from truth set (same-base filtering not needed)")
    same_base_negatives_df = pd.DataFrame()  # Empty - use explicit negatives
else:
    same_base_negatives_df = generate_same_base_negatives(tool_dfs, truth_pos, tool_base_cols)

# Always write the negative sites file (even if empty) to satisfy Snakemake output tracking
if not same_base_negatives_df.empty:
    print(f"\nTotal same-base negative sites: {len(same_base_negatives_df)}")
else:
    print("No same-base negative sites generated (using explicit negatives or empty)")
    # Create empty DataFrame with expected columns
    same_base_negatives_df = pd.DataFrame(columns=[
        'transcript', 'position', 'modification_type', 'ref_base', 'site_type'
    ])

same_base_negatives_df.to_csv(output_negatives, sep='\t', index=False)
print(f"Saved same-base negative sites to {output_negatives}")


# ── Compute Metrics with Same-Base Negatives ──────────────────────────────────

records = []
all_mod_types = set(truth_pos['modification_type'].unique())

for tool, pred_df in tool_dfs.items():
    for mod_type in all_mod_types:
        truth_subset = truth_pos[truth_pos['modification_type'] == mod_type].copy()

        if truth_subset.empty:
            continue

        same_base_neg_subset = same_base_negatives_df[
            same_base_negatives_df['modification_type'] == mod_type
        ].copy() if not same_base_negatives_df.empty else pd.DataFrame()

        all_tx = set(truth_subset['transcript'].unique())
        if not same_base_neg_subset.empty:
            all_tx.update(same_base_neg_subset['transcript'].unique())

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
            total_negative = 0
            if not same_base_neg_subset.empty:
                total_negative = len(same_base_neg_subset)
                for _, neg_row in same_base_neg_subset.iterrows():
                    tx = neg_row['transcript']
                    pos = neg_row['position']

                    # Check if there's any prediction within window
                    nearby_pred = pred_subset[
                        (pred_subset['transcript'] == tx) &
                        (pred_subset['position'].between(pos - window, pos + window))
                    ]

                    if nearby_pred.empty:
                        tn += 1

            called_sites = len(pred_subset)
            total_truth = len(truth_subset)

            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = (2 * precision * recall / (precision + recall)
                  if (precision + recall) > 0 else 0)

            # Specificity: TN / (TN + FP) or TN / total_negative
            if total_negative == 0:
                specificity = np.nan
            else:
                specificity = tn / total_negative

            # Matthews Correlation Coefficient
            if total_negative == 0:
                mcc = np.nan
            else:
                mcc = compute_mcc(tp, fp, fn, tn)

            # AUPRC and AUROC (only if score column exists)
            auprc, auroc = compute_ranking_metrics(
                pred_subset, truth_subset, same_base_neg_subset, score_col, window, tool_name=tool
            )

            records.append({
                "tool": tool,
                "modification_type": mod_type,
                "window": window,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "tn": tn,
                "precision": precision,
                "recall": recall,
                "f1": f1,
                "specificity": specificity,
                "mcc": mcc,
                "auprc": auprc,
                "auroc": auroc,
                "called_sites": called_sites,
                "total_truth": total_truth,
                "total_negative": total_negative,
                "negative_type": "same_base"
            })

# ── Write Output ──────────────────────────────────────────────────────────────

if records:
    out_df = pd.DataFrame(records).sort_values(
        ["tool", "modification_type", "window", "f1"],
        ascending=[True, True, True, False]
    )
else:
    out_df = pd.DataFrame(columns=[
        "tool", "modification_type", "window", "tp", "fp", "fn", "tn",
        "precision", "recall", "f1", "specificity", "mcc", "auprc", "auroc",
        "called_sites", "total_truth", "total_negative", "negative_type"
    ])

out_df.to_csv(output_metrics, sep='\t', index=False)
print(f"\nSaved same-base negative metrics to {output_metrics}")
print(f"Total records: {len(out_df)}")
