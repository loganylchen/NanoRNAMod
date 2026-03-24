"""
Per-tool internal score column optimization for RNA modification detection.

Features:
- Evaluates ALL available score columns for each tool
- For each score column, finds optimal threshold maximizing F1 score
- Computes AUROC and AUPRC for each score column
- Handles both p-value (lower=better) and probability (higher=better) scores
- Selects the best score column per tool based on F1 score

Output files:
- optimal_score_per_tool.tsv - best score column per tool with metrics
- all_scores_evaluation.tsv - detailed metrics for all score columns evaluated
"""

import os
import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_curve, roc_curve, auc


def tool_from_path(path):
    """Extract tool name from path."""
    parts = path.split(os.sep)
    for i, part in enumerate(parts):
        if part == "modifications" and i + 1 < len(parts):
            return parts[i + 1]
    return os.path.basename(os.path.dirname(os.path.dirname(path)))


def normalize_columns(df):
    """Map various column names to standard 'transcript' and 'position'."""
    col_mapping = {}

    transcript_cols = ['transcript', 'id', 'ref_id', 'chrom']
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


# Tool-specific score columns to evaluate (in priority order)
# Based on documented tool outputs from tools_comparison.md and format.py
# Priority: p-value/adjusted p-value first, then other metrics
TOOL_SCORE_COLUMNS = {
    # Signal-based comparative tools
    # xPore output: transcript, position, kmer, statistic, p_value, adjusted_p_value, direction
    'xpore': ['p_value', 'adjusted_p_value', 'statistic', 'pval', 'padj'],
    # Nanocompore output: transcript, position, ref_kmer, coverage, p_value, p_value_glm, p_value_ks, pass_glm, pass_ks
    'nanocompore': ['p_value_glm', 'p_value_ks', 'p_value', 'GMM_logit_pvalue', 'KS_dwell_pvalue', 'KS_intensity_pvalue'],
    # Baleen output: transcript, position, pvalue, padj, effect_size, mod_type
    'baleen': ['pvalue', 'padj', 'effect_size', 'score'],
    # pyBaleen output: transcript, position, pvalue, padj, effect_size
    'pybaleen': ['pvalue', 'padj', 'effect_size', 'p_value', 'adjusted_p_value', 'stoichiometry'],

    # Error-based comparative tools
    # DiffErr output (from format.py): -log10 P value, -log10 FDR, odds ratio, G statistic
    'differr': ['-log10 P value', '-log10 FDR', 'odds ratio', 'G statistic', 'score'],
    # DRUMMER output: various statistical columns
    'drummer': ['p_value', 'pvalue', 'padj', 'OR_padj', 'G_padj', 'odds_ratio', 'log2_(OR)'],
    # ELIGOS2 output: chr, start, end, ref_base, coverage, error_rate, p_value, q_value
    'eligos2': ['p_value', 'q_value', 'error_rate', 'pvalue', 'qvalue'],
    # EpiNano output (from format.py): chrom, pos, ref, strand, ko_feature, wt_feature, delta_sum_err, z_scores, z_score_prediction
    'epinano': ['z_scores', 'delta_sum_err', 'z_score_prediction', 'score'],
    # psipore output: transcript, position, psi_score, p_value, coverage
    'psipore': ['p_value', 'psi_score', 'pvalue', 'score'],

    # Single-sample ML/DL tools
    # TandemMod output: transcript, position, modification, probability, stoichiometry
    'tandemmod': ['probability', 'stoichiometry', 'prob', 'score'],
    # DirectRM output: transcript, position, modification, probability
    'directrm': ['probability', 'prob', 'score'],
    # m6ATM output: transcript, position, probability, stoichiometry
    'm6atm': ['probability', 'stoichiometry', 'prob', 'score'],
    # Rnano output: transcript, position, modification, probability
    'rnano': ['probability', 'prob', 'score'],
    # NanoPSU output: transcript, position, psi_probability, coverage
    'nanopsu': ['psi_probability', 'probability', 'pvalue', 'score'],
    # NanoMUD output: transcript, position, modification, probability, coverage
    'nanomud': ['probability', 'pvalue', 'p_value', 'score'],
    # Penguin output: transcript, position, psi_probability, score
    'penguin': ['psi_probability', 'probability', 'score', 'pvalue'],
}


def get_available_score_columns(df, tool_name=None):
    """
    Get all available score columns for a tool.

    Returns a list of (column_name, is_pvalue) tuples.
    Uses substring matching to handle columns with prefixes/suffixes.
    """
    available = []
    df_columns = list(df.columns)
    print(f"[DEBUG] get_available_score_columns for {tool_name}: columns = {df_columns}")

    if tool_name and tool_name in TOOL_SCORE_COLUMNS:
        potential_patterns = TOOL_SCORE_COLUMNS[tool_name]
        print(f"[DEBUG] {tool_name}: using tool-specific patterns: {potential_patterns}")
    else:
        # Generic score columns if tool not recognized
        potential_patterns = [
            'p_value', 'pvalue', 'pval', '-log10_pvalue', '-log10 P value',
            'score', 'probability', 'prob', 'mod_score', 'mod_prob',
            'logit_pvalue', 'logit', 'z_score', 'z_scores', 'esb', 'oddsR',
            'diff_mod', 'diff_mod_frac', 'mod_ratio', 'stoichiometry',
            'kmer_score', 'delta_sum_err', 'z_score_prediction', 'coverage',
            'psi_probability', 'psi_score', 'effect_size', 'mean_p_mod',
            'adjusted_p_value', 'padj', 'q_value', 'statistic',
        ]

    def match_column(pattern):
        """Match a pattern against dataframe columns using multiple strategies."""
        pattern_lower = pattern.lower().replace(' ', '_')

        # 1. Exact match
        if pattern in df_columns:
            return pattern

        # 2. Case-insensitive exact match
        for col in df_columns:
            if col.lower().replace(' ', '_') == pattern_lower:
                return col

        # 3. Substring match (pattern is contained in column name)
        for col in df_columns:
            col_lower = col.lower().replace(' ', '_')
            if pattern_lower in col_lower:
                return col

        # 4. Reverse substring match (column name is contained in pattern)
        for col in df_columns:
            col_lower = col.lower().replace(' ', '_')
            if col_lower in pattern_lower:
                return col

        return None

    for pattern in potential_patterns:
        matched = match_column(pattern)
        if matched and matched not in [a[0] for a in available]:
            is_pval = is_pvalue_column(matched, tool_name)
            available.append((matched, is_pval))

    return available


def is_pvalue_column(col_name, tool_name=None):
    """
    Determine if a column is a p-value (lower=better) or score (higher=better).

    Returns:
        True if lower values indicate more significant findings (p-values)
        False if higher values indicate more significant findings (scores, probabilities)
    """
    if not col_name:
        return False

    col_lower = col_name.lower().replace(' ', '_')

    # -log10 columns are higher-is-better (negative log of p-value)
    if '-log10' in col_lower:
        return False

    # Explicit p-value indicators (lower is better)
    if any(x in col_lower for x in ['p_value', 'pvalue', 'pval', 'padj', 'adjusted_p_value', 'fdr', 'q_value', 'qval']):
        return True

    # Probability columns (higher = more modification, better)
    if any(x in col_lower for x in ['probability', 'prob', 'psi_probability', 'mean_p_mod']):
        return False

    # Stoichiometry (higher = more modification)
    if 'stoichiometry' in col_lower:
        return False

    # Effect size (higher = larger difference = more modification)
    if 'effect_size' in col_lower:
        return False

    # Logit values: higher means more significant
    if 'logit' in col_lower:
        return False

    # Z-scores: absolute higher means more significant
    if 'z_score' in col_lower:
        return False

    # Score columns typically higher=better
    if 'score' in col_lower or 'psi_score' in col_lower:
        return False

    # Error-based metrics (delta, diff): higher absolute = more difference
    if 'delta' in col_lower or 'diff_mod' in col_lower or 'mod_ratio' in col_lower:
        return False

    # Default: assume higher is better
    return False


def generate_thresholds(scores, n_thresholds=100):
    """
    Generate threshold values based on score distribution.
    """
    scores = scores.dropna()

    if scores.empty:
        return []

    # Generate percentile-based thresholds for even coverage
    percentiles = np.linspace(0, 100, n_thresholds)
    percentile_thresholds = np.percentile(scores, percentiles)

    # Get unique values
    thresholds = sorted(set(percentile_thresholds))

    return thresholds


def compute_metrics_at_threshold(pred_df, truth_pos, threshold, score_col,
                                 is_pvalue, window=0):
    """
    Compute precision, recall, F1 at a given threshold.
    """
    # Filter predictions by threshold
    if is_pvalue:
        filtered = pred_df[pred_df[score_col] <= threshold].copy()
    else:
        filtered = pred_df[pred_df[score_col] >= threshold].copy()

    if filtered.empty:
        return {
            'threshold': threshold,
            'tp': 0, 'fp': 0, 'fn': len(truth_pos),
            'precision': 0.0, 'recall': 0.0, 'f1': 0.0,
            'predictions': 0
        }

    # Count TP (predictions within window of truth sites)
    matched_pred_indices = set()
    tp = 0

    for _, truth_row in truth_pos.iterrows():
        tx = truth_row['transcript']
        pos = truth_row['position']

        matches = filtered[
            (filtered['transcript'] == tx) &
            (filtered['position'].between(pos - window, pos + window))
        ]

        if not matches.empty:
            tp += 1
            matched_pred_indices.update(matches.index.tolist())

    fp = len(filtered) - len(matched_pred_indices)
    fn = len(truth_pos) - tp

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)

    return {
        'threshold': threshold,
        'tp': tp, 'fp': fp, 'fn': fn,
        'precision': precision, 'recall': recall, 'f1': f1,
        'predictions': len(filtered)
    }


def compute_auc_metrics(pred_df, truth_pos, score_col, is_pvalue, window=0):
    """
    Compute AUROC and AUPRC for a score column.

    Creates binary labels based on whether predictions are within window
    of truth sites.
    """
    if pred_df.empty or score_col not in pred_df.columns:
        return {'auprc': np.nan, 'auroc': np.nan}

    scores = pred_df[score_col].dropna()
    if scores.empty:
        return {'auprc': np.nan, 'auroc': np.nan}

    # Create binary labels: 1 if within window of a truth site, else 0
    labels = np.zeros(len(pred_df))

    # Create a mapping for truth sites for faster lookup
    truth_dict = {}
    for _, row in truth_pos.iterrows():
        tx = row['transcript']
        pos = row['position']
        if tx not in truth_dict:
            truth_dict[tx] = []
        truth_dict[tx].append(pos)

    for idx, row in pred_df.iterrows():
        tx = row['transcript']
        pos = row['position']

        if tx in truth_dict:
            for truth_pos in truth_dict[tx]:
                if abs(pos - truth_pos) <= window:
                    labels[idx] = 1
                    break

    # Filter to rows with valid scores
    valid_mask = pred_df[score_col].notna()
    y_true = labels[valid_mask]
    y_scores = pred_df.loc[valid_mask, score_col].values

    if len(np.unique(y_true)) < 2:
        # Need both positive and negative samples for ROC/PR
        return {'auprc': np.nan, 'auroc': np.nan}

    # For p-values (lower=better), invert scores for AUROC/AUPRC computation
    # sklearn expects higher scores = more likely positive
    if is_pvalue:
        y_scores = -y_scores

    try:
        # Compute AUPRC
        precision, recall, _ = precision_recall_curve(y_true, y_scores)
        auprc = auc(recall, precision)

        # Compute AUROC
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        auroc = auc(fpr, tpr)

        return {'auprc': auprc, 'auroc': auroc}
    except Exception:
        return {'auprc': np.nan, 'auroc': np.nan}


def evaluate_score_column(pred_df, truth_pos, tool_name, mod_type,
                          score_col, is_pvalue, window=0, n_thresholds=100):
    """
    Evaluate a single score column across multiple thresholds.

    Returns:
        dict with optimal_threshold, optimal_f1, optimal_precision,
              optimal_recall, auprc, auroc, and all_evaluation_records
    """
    if pred_df.empty or score_col not in pred_df.columns:
        return {
            'score_column': score_col,
            'score_type': 'p-value' if is_pvalue else 'probability',
            'optimal_threshold': np.nan,
            'optimal_f1': 0.0,
            'optimal_precision': 0.0,
            'optimal_recall': 0.0,
            'auprc': np.nan,
            'auroc': np.nan,
            'evaluation_records': []
        }

    # Generate thresholds
    scores = pred_df[score_col]
    thresholds = generate_thresholds(scores, n_thresholds=n_thresholds)

    if not thresholds:
        thresholds = [scores.median()]

    # Compute AUC metrics
    auc_metrics = compute_auc_metrics(pred_df, truth_pos, score_col, is_pvalue, window)

    # Evaluate at each threshold
    best_f1 = -1
    best_record = None
    eval_records = []

    for thresh in thresholds:
        metrics = compute_metrics_at_threshold(
            pred_df, truth_pos, thresh, score_col, is_pvalue, window
        )

        record = {
            'tool': tool_name,
            'modification_type': mod_type,
            'score_column': score_col,
            'score_type': 'p-value' if is_pvalue else 'probability',
            'threshold': thresh,
            **metrics
        }
        eval_records.append(record)

        if metrics['f1'] >= best_f1:
            best_f1 = metrics['f1']
            best_record = record

    if best_record is None:
        best_record = {
            'optimal_threshold': np.nan,
            'optimal_f1': 0.0,
            'optimal_precision': 0.0,
            'optimal_recall': 0.0,
        }
    else:
        best_record = {
            'optimal_threshold': best_record['threshold'],
            'optimal_f1': best_record['f1'],
            'optimal_precision': best_record['precision'],
            'optimal_recall': best_record['recall'],
        }

    return {
        'score_column': score_col,
        'score_type': 'p-value' if is_pvalue else 'probability',
        **best_record,
        **auc_metrics,
        'evaluation_records': eval_records
    }


def evaluate_all_score_columns(pred_df, truth_pos, tool_name, mod_type,
                               window=0, n_thresholds=100):
    """
    Evaluate all available score columns for a tool.

    Returns:
        list of evaluation results for each score column
    """
    available_columns = get_available_score_columns(pred_df, tool_name)

    if not available_columns:
        print(f"Warning: No score columns found for tool {tool_name}")
        return []

    results = []

    for score_col, is_pvalue in available_columns:
        result = evaluate_score_column(
            pred_df, truth_pos, tool_name, mod_type,
            score_col, is_pvalue, window, n_thresholds
        )
        results.append(result)

    return results


def main():
    # When called by Snakemake
    if 'snakemake' in globals():
        result_files = snakemake.input.results
        truth_set_path = snakemake.input.truth_set
        output_optimal = snakemake.output.optimal
        output_all = snakemake.output.all_eval

        window_param = snakemake.params.get('window', 0)
        n_thresholds = snakemake.params.get('n_thresholds', 100)

        if isinstance(window_param, list):
            windows = window_param
        else:
            windows = [int(window_param)]
    else:
        import argparse
        parser = argparse.ArgumentParser(
            description='Score column optimization for RNA modification detection'
        )
        parser.add_argument('--results', nargs='+', required=True,
                          help='Tool result TSV files')
        parser.add_argument('--truth', required=True,
                          help='Ground truth TSV file')
        parser.add_argument('--output-optimal', required=True,
                          help='Output file for optimal score per tool')
        parser.add_argument('--output-all', required=True,
                          help='Output file for all score columns evaluated')
        parser.add_argument('--window', type=int, default=0,
                          help='Positional tolerance window')
        parser.add_argument('--n-thresholds', type=int, default=100,
                          help='Number of thresholds to evaluate per score column')
        args = parser.parse_args()

        result_files = args.results
        truth_set_path = args.truth
        output_optimal = args.output_optimal
        output_all = args.output_all
        windows = [args.window]
        n_thresholds = args.n_thresholds

    # Load truth set
    truth = pd.read_csv(truth_set_path, sep='\t')

    if 'label' in truth.columns:
        truth_pos = truth[truth['label'] != '-'].copy()
    else:
        truth_pos = truth.copy()

    if truth_pos.empty:
        print("Warning: No positive sites in truth set")
        # Write empty outputs
        pd.DataFrame(columns=[
            'tool', 'score_column', 'score_type', 'optimal_threshold',
            'optimal_f1', 'optimal_precision', 'optimal_recall',
            'auprc', 'auroc', 'all_available_scores'
        ]).to_csv(output_optimal, sep='\t', index=False)

        pd.DataFrame(columns=[
            'tool', 'modification_type', 'score_column', 'score_type',
            'threshold', 'tp', 'fp', 'fn', 'precision', 'recall', 'f1', 'predictions'
        ]).to_csv(output_all, sep='\t', index=False)
        return

    # Normalize truth columns
    truth_pos = normalize_columns(truth_pos)

    # Load tool results
    tool_dfs = {}

    for f in result_files:
        tool = tool_from_path(f)
        try:
            df = pd.read_csv(f, sep='\t')
            df = normalize_columns(df)

            if 'transcript' in df.columns and 'position' in df.columns:
                if tool not in tool_dfs:
                    tool_dfs[tool] = []
                tool_dfs[tool].append(df)
        except Exception as e:
            print(f"Warning: could not read {f}: {e}")

    # Concatenate results per tool
    for tool in tool_dfs:
        tool_dfs[tool] = pd.concat(tool_dfs[tool], ignore_index=True).drop_duplicates(
            subset=['transcript', 'position']
        )

    # Get modification types
    if 'modification_type' in truth_pos.columns:
        all_mod_types = truth_pos['modification_type'].unique()
    else:
        all_mod_types = ['unknown']

    # Evaluate all score columns for each tool
    all_evaluation_records = []
    optimal_per_tool = []

    for tool, pred_df in tool_dfs.items():
        # Get all available score columns for this tool
        available_cols = get_available_score_columns(pred_df, tool)
        all_score_names = [col for col, _ in available_cols]

        if not all_score_names:
            print(f"Warning: No score columns found for tool {tool}")
            continue

        for mod_type in all_mod_types:
            # Filter truth and predictions for this mod type
            if 'modification_type' in truth_pos.columns:
                truth_subset = truth_pos[
                    truth_pos['modification_type'] == mod_type
                ].copy()
            else:
                truth_subset = truth_pos.copy()

            if truth_subset.empty:
                continue

            all_tx = set(truth_subset['transcript'].unique())
            pred_subset = pred_df[pred_df['transcript'].isin(all_tx)].copy()

            if pred_subset.empty:
                continue

            for window in windows:
                # Evaluate all score columns
                score_results = evaluate_all_score_columns(
                    pred_subset, truth_subset, tool, mod_type,
                    window, n_thresholds
                )

                # Add all evaluation records
                for result in score_results:
                    all_evaluation_records.extend(result.get('evaluation_records', []))

                # Find the best score column for this tool/mod_type
                if score_results:
                    # Sort by F1, then AUROC
                    scored_results = [
                        (r, r.get('optimal_f1', 0), r.get('auprc', 0))
                        for r in score_results
                    ]
                    scored_results.sort(key=lambda x: (x[1], x[2]), reverse=True)

                    best = scored_results[0][0]

                    optimal_record = {
                        'tool': tool,
                        'modification_type': mod_type,
                        'window': window,
                        'score_column': best['score_column'],
                        'score_type': best['score_type'],
                        'optimal_threshold': best['optimal_threshold'],
                        'optimal_f1': best['optimal_f1'],
                        'optimal_precision': best['optimal_precision'],
                        'optimal_recall': best['optimal_recall'],
                        'auprc': best['auprc'],
                        'auroc': best['auroc'],
                        'all_available_scores': ','.join(all_score_names)
                    }
                    optimal_per_tool.append(optimal_record)

    # Write outputs
    if all_evaluation_records:
        all_df = pd.DataFrame(all_evaluation_records)
        all_df = all_df.sort_values(
            ['tool', 'modification_type', 'score_column', 'f1'],
            ascending=[True, True, True, False]
        )
    else:
        all_df = pd.DataFrame(columns=[
            'tool', 'modification_type', 'score_column', 'score_type',
            'threshold', 'tp', 'fp', 'fn', 'precision', 'recall', 'f1', 'predictions'
        ])

    if optimal_per_tool:
        optimal_df = pd.DataFrame(optimal_per_tool)
        optimal_df = optimal_df.sort_values(
            ['modification_type', 'window', 'optimal_f1'],
            ascending=[True, True, False]
        )
    else:
        optimal_df = pd.DataFrame(columns=[
            'tool', 'modification_type', 'window', 'score_column', 'score_type',
            'optimal_threshold', 'optimal_f1', 'optimal_precision',
            'optimal_recall', 'auprc', 'auroc', 'all_available_scores'
        ])

    all_df.to_csv(output_all, sep='\t', index=False)
    optimal_df.to_csv(output_optimal, sep='\t', index=False)

    print(f"Written all score columns evaluation to {output_all}")
    print(f"Written optimal score per tool to {output_optimal}")

    # Print summary
    if not optimal_df.empty:
        print("\n=== Optimal Score Column Per Tool ===")
        for _, row in optimal_df.iterrows():
            print(f"{row['tool']}: {row['score_column']} "
                  f"(F1={row['optimal_f1']:.3f}, AUROC={row['auroc']:.3f})")


if __name__ == '__main__':
    main()
