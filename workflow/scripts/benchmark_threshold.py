"""
Optimal threshold analysis for RNA modification detection tools.

Finds the optimal score threshold for each tool that maximizes F1 score
or balances precision/recall based on specified criteria.
"""

import os
import pandas as pd
import numpy as np


def tool_from_path(path):
    """Extract tool name from path like .../modifications/{tool}/{comp}/{tool}_results.tsv"""
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


def detect_score_column(df, tool_name=None):
    """Detect the most likely score column for ranking predictions."""
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
        'pybaleen': ['mod_score', 'score'],
    }

    if tool_name and tool_name in tool_score_map:
        for col in tool_score_map[tool_name]:
            if col in df.columns:
                return col
            for df_col in df.columns:
                if df_col.lower() == col.lower():
                    return df_col

    generic_score_cols = [
        'p_value', 'pvalue', 'pval', '-log10_pvalue', '-log10 P value',
        'score', 'probability', 'prob', 'mod_score', 'mod_prob',
        'logit_pvalue', 'logit', 'z_score', 'z_scores'
    ]

    for col in generic_score_cols:
        if col in df.columns:
            return col
        for df_col in df.columns:
            if df_col.lower().replace(' ', '_') == col.lower().replace(' ', '_'):
                return df_col

    return None


def compute_metrics_at_threshold(pred_df, truth_pos, threshold, score_col,
                                 is_pvalue, window=0):
    """
    Compute precision, recall, F1 at a given threshold.

    Args:
        pred_df: Predictions dataframe
        truth_pos: Positive truth sites (with transcript, position)
        threshold: Score threshold for calling positive
        score_col: Column name containing scores
        is_pvalue: Whether score_col is a p-value (lower=better) or score (higher=better)
        window: Positional tolerance window

    Returns:
        dict with tp, fp, fn, precision, recall, f1
    """
    # Filter predictions by threshold
    if is_pvalue:
        filtered = pred_df[pred_df[score_col] <= threshold].copy()
    else:
        filtered = pred_df[pred_df[score_col] >= threshold].copy()

    if filtered.empty:
        return {'tp': 0, 'fp': 0, 'fn': len(truth_pos),
                'precision': 0.0, 'recall': 0.0, 'f1': 0.0}

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

    return {'tp': tp, 'fp': fp, 'fn': fn,
            'precision': precision, 'recall': recall, 'f1': f1}


def find_optimal_threshold(pred_df, truth_pos, score_col, is_pvalue,
                          window=0, criterion='f1', n_thresholds=100):
    """
    Find optimal threshold by evaluating metrics at multiple thresholds.

    Args:
        pred_df: Predictions with score column
        truth_pos: Positive truth sites
        score_col: Score column name
        is_pvalue: Whether lower score is better
        window: Positional tolerance
        criterion: 'f1', 'precision', 'recall', 'balanced'
        n_thresholds: Number of thresholds to test

    Returns:
        dict with optimal_threshold, optimal_f1, and metrics at threshold
    """
    if pred_df.empty or score_col not in pred_df.columns:
        return None

    # Get valid scores
    scores = pred_df[score_col].dropna()

    if scores.empty:
        return None

    # Generate thresholds
    score_min, score_max = scores.min(), scores.max()
    thresholds = np.linspace(score_min, score_max, n_thresholds)

    best_result = None
    best_value = -np.inf

    for thresh in thresholds:
        metrics = compute_metrics_at_threshold(
            pred_df, truth_pos, thresh, score_col, is_pvalue, window
        )

        # Select value to optimize
        if criterion == 'f1':
            value = metrics['f1']
        elif criterion == 'precision':
            value = metrics['precision']
        elif criterion == 'recall':
            value = metrics['recall']
        elif criterion == 'balanced':
            # Base F1 score with balance penalty applied below
            value = metrics['f1']
        else:
            value = metrics['f1']

        # For balanced optimization, also consider precision-recall balance
        if criterion == 'balanced':
            # Prefer balanced precision/recall over pure F1
            balance_penalty = abs(metrics['precision'] - metrics['recall'])
            value = value - 0.1 * balance_penalty

        if value > best_value:
            best_value = value
            best_result = {
                'threshold': thresh,
                'f1': metrics['f1'],
                'precision': metrics['precision'],
                'recall': metrics['recall'],
                'tp': metrics['tp'],
                'fp': metrics['fp'],
                'fn': metrics['fn'],
            }

    return best_result


def main():
    # When called by Snakemake
    if 'snakemake' in globals():
        result_files = snakemake.input.results
        truth_set_path = snakemake.input.truth_set
        output_file = snakemake.output[0]
        window_param = snakemake.params.get('window', 0)
        criterion = snakemake.params.get('criterion', 'f1')

        if isinstance(window_param, list):
            windows = window_param
        else:
            windows = [int(window_param)]
    else:
        # Command line interface
        import argparse
        import sys

        parser = argparse.ArgumentParser(
            description='Find optimal score thresholds for modification detection tools'
        )
        parser.add_argument('--results', nargs='+', required=True,
                          help='Tool result TSV files')
        parser.add_argument('--truth', required=True,
                          help='Ground truth TSV file')
        parser.add_argument('--output', required=True,
                          help='Output TSV file')
        parser.add_argument('--window', type=int, default=0,
                          help='Positional tolerance window')
        parser.add_argument('--criterion', default='f1',
                          choices=['f1', 'precision', 'recall', 'balanced'],
                          help='Optimization criterion')
        args = parser.parse_args()

        result_files = args.results
        truth_set_path = args.truth
        output_file = args.output
        windows = [args.window]
        criterion = args.criterion

    # Load truth set
    truth = pd.read_csv(truth_set_path, sep='\t')

    if 'label' in truth.columns:
        truth_pos = truth[truth['label'] != '-'].copy()
    else:
        truth_pos = truth.copy()

    if truth_pos.empty:
        print("Warning: No positive sites in truth set")
        # Write empty output
        pd.DataFrame(columns=[
            'tool', 'modification_type', 'window', 'criterion',
            'optimal_threshold', 'optimal_f1', 'precision', 'recall',
            'tp', 'fp', 'fn'
        ]).to_csv(output_file, sep='\t', index=False)
        return

    # Load tool results
    tool_dfs = {}
    tool_score_cols = {}

    for f in result_files:
        tool = tool_from_path(f)
        try:
            df = pd.read_csv(f, sep='\t')
            df = normalize_columns(df)

            if 'transcript' in df.columns and 'position' in df.columns:
                score_col = detect_score_column(df, tool)
                if score_col:
                    if tool not in tool_score_cols:
                        tool_score_cols[tool] = score_col
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

    # Find optimal thresholds
    records = []
    all_mod_types = truth_pos['modification_type'].unique() if 'modification_type' in truth_pos.columns else ['unknown']

    for tool, pred_df in tool_dfs.items():
        score_col = tool_score_cols.get(tool)

        if not score_col or score_col not in pred_df.columns:
            print(f"Warning: No score column found for {tool}")
            continue

        # Check if p-value column
        score_col_lower = score_col.lower().replace(' ', '_')
        is_pvalue = any(x in score_col_lower for x in ['p_value', 'pvalue', 'pval'])

        for mod_type in all_mod_types:
            truth_subset = truth_pos[
                truth_pos['modification_type'] == mod_type
            ].copy() if 'modification_type' in truth_pos.columns else truth_pos

            if truth_subset.empty:
                continue

            # Get unique transcripts for this mod type
            all_tx = set(truth_subset['transcript'].unique())
            pred_subset = pred_df[pred_df['transcript'].isin(all_tx)].copy()

            if pred_subset.empty:
                continue

            for window in windows:
                result = find_optimal_threshold(
                    pred_subset, truth_subset, score_col, is_pvalue,
                    window, criterion
                )

                if result:
                    records.append({
                        'tool': tool,
                        'modification_type': mod_type,
                        'window': window,
                        'criterion': criterion,
                        'optimal_threshold': result['threshold'],
                        'optimal_f1': result['f1'],
                        'precision': result['precision'],
                        'recall': result['recall'],
                        'tp': result['tp'],
                        'fp': result['fp'],
                        'fn': result['fn'],
                    })

    # Write output
    if records:
        out_df = pd.DataFrame(records).sort_values(
            ['modification_type', 'window', 'optimal_f1'],
            ascending=[True, True, False]
        )
    else:
        out_df = pd.DataFrame(columns=[
            'tool', 'modification_type', 'window', 'criterion',
            'optimal_threshold', 'optimal_f1', 'precision', 'recall',
            'tp', 'fp', 'fn'
        ])

    out_df.to_csv(output_file, sep='\t', index=False)
    print(f"Written optimal thresholds to {output_file}")


if __name__ == '__main__':
    main()
