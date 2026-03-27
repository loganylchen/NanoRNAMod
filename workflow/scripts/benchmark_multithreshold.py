"""
Multi-threshold evaluation for RNA modification detection benchmarking.

Features:
- Automatic threshold generation from score distribution percentiles
- Custom threshold support
- Handles both p-value (lower=better) and probability (higher=better) scores
- Computes P/R/F1 at each threshold
- Finds optimal threshold per tool/modification type
- Cross-tool comparison using optimal thresholds

Output files:
- threshold_evaluation.tsv - metrics at each threshold
- optimal_thresholds.tsv - best threshold per tool/modification
- score_distributions.tsv - score statistics per tool
"""

import os
import pandas as pd
import numpy as np

# Import shared utilities
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_utils import tool_from_path, comparison_from_path, normalize_columns


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


def is_pvalue_column(col_name, tool_name=None):
    """
    Determine if a column contains raw p-values that need -log10 transformation.

    Returns True if the column is a raw p-value (smaller=better) that should
    be transformed to -log10(p-value) for uniform "higher=better" comparison.

    Returns False if the column is already:
    - A probability/score (higher=better)
    - Already -log10 transformed (higher=better)
    """
    if not col_name:
        return False

    col_lower = col_name.lower().replace(' ', '_')

    # Already -log10 transformed columns are higher=better, no transformation needed
    if '-log10' in col_lower:
        return False

    # Explicit p-value indicators (raw p-values, need transformation)
    if any(x in col_lower for x in ['p_value', 'pvalue', 'pval', 'fdr', 'qval', 'padj', 'adj_p']):
        return True

    # Tool-specific heuristics
    if tool_name == 'xpore' and 'p_value' in col_lower:
        return True
    if tool_name == 'nanocompore' and 'pvalue' in col_lower:
        return True
    if tool_name == 'differr' and '-log10' in col_lower:
        return False  # -log10 is already higher=better
    if tool_name in ['tandemmod', 'directrm', 'm6atm', 'rnano'] and 'probability' in col_lower:
        return False

    # Default: assume higher is better for score/probability columns
    return False


def transform_pvalue_to_log10(scores, is_pvalue):
    """
    Transform p-values to -log10(p-values) for uniform "higher=better" comparison.

    Args:
        scores: Series of score values
        is_pvalue: Whether these are raw p-values

    Returns:
        Transformed scores (if is_pvalue) or original scores
    """
    if not is_pvalue:
        return scores

    # Convert to numeric
    scores_numeric = pd.to_numeric(scores, errors='coerce')

    # Transform: -log10(p-value), handling edge cases
    # Clamp p-values to avoid log(0)
    scores_numeric = scores_numeric.clip(lower=1e-300)

    return -np.log10(scores_numeric)


def generate_thresholds(scores, n_thresholds=50, min_thresholds=10, custom_thresholds=None):
    """
    Generate threshold values based on score distribution.

    Uses percentiles to ensure even coverage across the score range.
    """
    # Convert to numeric, coercing errors to NaN, then drop NAs
    scores = pd.to_numeric(scores, errors='coerce').dropna()

    if scores.empty:
        return []

    # Generate percentile-based thresholds
    percentiles = np.linspace(0, 100, n_thresholds)
    percentile_thresholds = np.percentile(scores, percentiles)

    # Add unique values only
    thresholds = sorted(set(percentile_thresholds))

    # Ensure minimum number of thresholds
    if len(thresholds) < min_thresholds:
        # Create evenly spaced thresholds between min and max
        min_val, max_val = scores.min(), scores.max()
        thresholds = np.linspace(min_val, max_val, min_thresholds).tolist()

    # Add custom thresholds if provided
    if custom_thresholds:
        thresholds.extend([t for t in custom_thresholds
                         if t >= scores.min() and t <= scores.max()])

    # Remove duplicates and sort
    thresholds = sorted(set(thresholds))

    return thresholds


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
    # Convert score column to numeric for comparison
    pred_df = pred_df.copy()
    pred_df['_numeric_score'] = pd.to_numeric(pred_df[score_col], errors='coerce')

    # Filter predictions by threshold
    if is_pvalue:
        filtered = pred_df[pred_df['_numeric_score'] <= threshold].copy()
    else:
        filtered = pred_df[pred_df['_numeric_score'] >= threshold].copy()

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


def compute_score_distribution(pred_df, score_col):
    """Compute summary statistics for score distribution."""
    # Convert to numeric, coercing errors to NaN, then drop NAs
    scores = pd.to_numeric(pred_df[score_col], errors='coerce').dropna()

    if scores.empty:
        return {
            'count': 0,
            'min': np.nan, 'max': np.nan, 'mean': np.nan, 'std': np.nan,
            'median': np.nan, 'q25': np.nan, 'q75': np.nan,
        }

    return {
        'count': len(scores),
        'min': scores.min(),
        'max': scores.max(),
        'mean': scores.mean(),
        'std': scores.std(),
        'median': scores.median(),
        'q25': scores.quantile(0.25),
        'q75': scores.quantile(0.75),
    }


def evaluate_tool_at_thresholds(pred_df, truth_pos, tool_name, mod_type,
                                score_col, is_pvalue, thresholds, window=0):
    """
    Evaluate a tool at multiple thresholds.

    Returns:
        tuple: (evaluation_records, optimal_record, score_dist_record)
    """
    if pred_df.empty or not score_col or score_col not in pred_df.columns:
        return [], None, None

    if not thresholds:
        return [], None, None

    records = []
    best_f1 = -1
    best_record = None

    for thresh in thresholds:
        metrics = compute_metrics_at_threshold(
            pred_df, truth_pos, thresh, score_col, is_pvalue, window
        )

        record = {
            'tool': tool_name,
            'modification_type': mod_type,
            'window': window,
            'score_column': score_col,
            'is_pvalue': is_pvalue,
            **metrics
        }
        records.append(record)

        # Track optimal threshold (first record with highest F1)
        if metrics['f1'] >= best_f1:
            best_f1 = metrics['f1']
            best_record = record.copy()
            best_record['criterion'] = 'f1'

    # Fallback: if no positive matches found, use median score threshold
    if best_f1 == 0 and records:
        median_score = pd.to_numeric(pred_df[score_col], errors='coerce').median()
        print(f"Warning: No positive F1 for {tool_name}/{mod_type}. Using median score: {median_score}")
        # Use the median as fallback threshold
        metrics = compute_metrics_at_threshold(
            pred_df, truth_pos, median_score, score_col, is_pvalue, window
        )
        best_record = {
            'tool': tool_name,
            'modification_type': mod_type,
            'window': window,
            'score_column': score_col,
            'is_pvalue': is_pvalue,
            **metrics
        }
        best_record['criterion'] = 'f1_median_fallback'

    # Score distribution
    score_dist = compute_score_distribution(pred_df, score_col)
    score_dist_record = {
        'tool': tool_name,
        'modification_type': mod_type,
        'score_column': score_col,
        'is_pvalue': is_pvalue,
        **score_dist
    }

    return records, best_record, score_dist_record


def main():
    # When called by Snakemake
    if 'snakemake' in globals():
        result_files = snakemake.input.results
        truth_set_path = snakemake.input.truth_set
        output_eval = snakemake.output.evaluation
        output_optimal = snakemake.output.optimal
        output_dist = snakemake.output.distribution

        window_param = snakemake.params.get('window', 0)
        n_thresholds = snakemake.params.get('n_thresholds', 50)
        custom_thresholds = snakemake.params.get('custom_thresholds', [])

        # Handle custom_thresholds as string (from YAML config)
        if isinstance(custom_thresholds, str):
            import ast
            try:
                custom_thresholds = ast.literal_eval(custom_thresholds)
            except (ValueError, SyntaxError):
                custom_thresholds = []

        if isinstance(window_param, list):
            windows = window_param
        else:
            windows = [int(window_param)]
    else:
        import argparse
        parser = argparse.ArgumentParser(
            description='Multi-threshold evaluation for benchmarking'
        )
        parser.add_argument('--results', nargs='+', required=True,
                          help='Tool result TSV files')
        parser.add_argument('--truth', required=True,
                          help='Ground truth TSV file')
        parser.add_argument('--output-eval', required=True,
                          help='Output file for threshold evaluation')
        parser.add_argument('--output-optimal', required=True,
                          help='Output file for optimal thresholds')
        parser.add_argument('--output-dist', required=True,
                          help='Output file for score distributions')
        parser.add_argument('--window', type=int, default=0,
                          help='Positional tolerance window')
        parser.add_argument('--n-thresholds', type=int, default=50,
                          help='Number of thresholds to evaluate')
        parser.add_argument('--custom-thresholds', nargs='+', type=float,
                          help='Custom thresholds to evaluate')
        args = parser.parse_args()

        result_files = args.results
        truth_set_path = args.truth
        output_eval = args.output_eval
        output_optimal = args.output_optimal
        output_dist = args.output_dist
        windows = [args.window]
        n_thresholds = args.n_thresholds
        custom_thresholds = args.custom_thresholds or []

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
            'tool', 'modification_type', 'window', 'score_column', 'is_pvalue',
            'threshold', 'tp', 'fp', 'fn', 'precision', 'recall', 'f1', 'predictions'
        ]).to_csv(output_eval, sep='\t', index=False)

        pd.DataFrame(columns=[
            'tool', 'modification_type', 'window', 'score_column', 'is_pvalue',
            'threshold', 'f1', 'precision', 'recall', 'tp', 'fp', 'fn', 'criterion'
        ]).to_csv(output_optimal, sep='\t', index=False)

        pd.DataFrame(columns=[
            'tool', 'modification_type', 'score_column', 'is_pvalue',
            'count', 'min', 'max', 'mean', 'std', 'median', 'q25', 'q75'
        ]).to_csv(output_dist, sep='\t', index=False)
        return

    # Load tool results
    tool_dfs = {}
    tool_score_cols = {}
    tool_is_pvalue = {}
    tool_comparisons = {}  # Track comparisons per tool

    for f in result_files:
        tool = tool_from_path(f)
        comparison = comparison_from_path(f)
        try:
            df = pd.read_csv(f, sep='\t')
            df = normalize_columns(df)

            # Add comparison column to track source
            df['_comparison'] = comparison

            if 'transcript' in df.columns and 'position' in df.columns:
                score_col = detect_score_column(df, tool)
                if score_col:
                    if tool not in tool_dfs:
                        tool_dfs[tool] = []
                        tool_score_cols[tool] = score_col
                        tool_is_pvalue[tool] = is_pvalue_column(score_col, tool)
                        tool_comparisons[tool] = set()
                    tool_dfs[tool].append(df)
                    tool_comparisons[tool].add(comparison)
        except Exception as e:
            print(f"Warning: could not read {f}: {e}")

    # Concatenate results per tool
    for tool in tool_dfs:
        tool_dfs[tool] = pd.concat(tool_dfs[tool], ignore_index=True).drop_duplicates(
            subset=['transcript', 'position']
        )

    # Evaluate at multiple thresholds
    all_eval_records = []
    all_optimal_records = []
    all_dist_records = []

    all_mod_types = truth_pos['modification_type'].unique() if 'modification_type' in truth_pos.columns else ['unknown']

    for tool, pred_df in tool_dfs.items():
        score_col = tool_score_cols.get(tool)
        is_pvalue = tool_is_pvalue.get(tool, False)

        if not score_col:
            continue

        for mod_type in all_mod_types:
            truth_subset = truth_pos[
                truth_pos['modification_type'] == mod_type
            ].copy() if 'modification_type' in truth_pos.columns else truth_pos

            if truth_subset.empty:
                continue

            all_tx = set(truth_subset['transcript'].unique())
            pred_subset = pred_df[pred_df['transcript'].isin(all_tx)].copy()

            if pred_subset.empty or score_col not in pred_subset.columns:
                continue

            # Transform p-values to -log10 scale for uniform "higher=better" comparison
            if is_pvalue:
                # Create a new column for transformed scores
                transformed_col = f'{score_col}_log10'
                # Convert to numeric and apply -log10 transformation
                pred_subset[transformed_col] = pd.to_numeric(
                    pred_subset[score_col], errors='coerce'
                ).clip(lower=1e-300)
                # Apply -log10 transformation
                pred_subset[transformed_col] = -np.log10(pred_subset[transformed_col])

                # Generate thresholds from transformed scores (higher is better)
                transformed_scores = pred_subset[transformed_col].dropna()
                thresholds = generate_thresholds(transformed_scores, n_thresholds=n_thresholds,
                                                custom_thresholds=custom_thresholds)

                # Use transformed column for evaluation (is_pvalue=False after transformation)
                for window in windows:
                    eval_records, optimal_record, dist_record = evaluate_tool_at_thresholds(
                        pred_subset, truth_subset, tool, mod_type,
                        transformed_col, False,  # is_pvalue=False after transformation
                        thresholds, window
                    )

                    # Add original score column info for reporting
                    for rec in eval_records:
                        rec['original_score_column'] = score_col
                        # Convert threshold back to original p-value for reference
                        rec['original_threshold'] = 10 ** (-rec['threshold']) if rec['threshold'] else np.nan

                    all_eval_records.extend(eval_records)
                    if optimal_record:
                        optimal_record['original_score_column'] = score_col
                        optimal_record['original_threshold'] = 10 ** (-optimal_record['threshold']) if optimal_record['threshold'] else np.nan
                        all_optimal_records.append(optimal_record)
                    if dist_record:
                        dist_record['window'] = window
                        dist_record['original_score_column'] = score_col
                        all_dist_records.append(dist_record)
            else:
                # Not a p-value column, use original scores
                thresholds = generate_thresholds(pred_subset[score_col], n_thresholds=n_thresholds,
                                                custom_thresholds=custom_thresholds)

                for window in windows:
                    eval_records, optimal_record, dist_record = evaluate_tool_at_thresholds(
                        pred_subset, truth_subset, tool, mod_type,
                        score_col, False, thresholds, window
                    )

                    all_eval_records.extend(eval_records)
                    if optimal_record:
                        all_optimal_records.append(optimal_record)
                    if dist_record:
                        dist_record['window'] = window
                        all_dist_records.append(dist_record)

    # Write outputs
    if all_eval_records:
        eval_df = pd.DataFrame(all_eval_records)
        # Sort by tool, mod_type, f1 descending
        eval_df = eval_df.sort_values(
            ['tool', 'modification_type', 'window', 'f1'],
            ascending=[True, True, True, False]
        )
    else:
        eval_df = pd.DataFrame(columns=[
            'tool', 'modification_type', 'window', 'score_column', 'is_pvalue',
            'threshold', 'tp', 'fp', 'fn', 'precision', 'recall', 'f1', 'predictions',
            'original_score_column', 'original_threshold'
        ])

    if all_optimal_records:
        optimal_df = pd.DataFrame(all_optimal_records)
        optimal_df = optimal_df.sort_values(
            ['modification_type', 'window', 'f1'],
            ascending=[True, True, False]
        )
    else:
        optimal_df = pd.DataFrame(columns=[
            'tool', 'modification_type', 'window', 'score_column', 'is_pvalue',
            'threshold', 'f1', 'precision', 'recall', 'tp', 'fp', 'fn', 'criterion',
            'original_score_column', 'original_threshold'
        ])

    if all_dist_records:
        dist_df = pd.DataFrame(all_dist_records)
        dist_df = dist_df.sort_values(['tool', 'modification_type'])
    else:
        dist_df = pd.DataFrame(columns=[
            'tool', 'modification_type', 'score_column', 'is_pvalue',
            'count', 'min', 'max', 'mean', 'std', 'median', 'q25', 'q75',
            'window', 'original_score_column'
        ])

    eval_df.to_csv(output_eval, sep='\t', index=False)
    optimal_df.to_csv(output_optimal, sep='\t', index=False)
    dist_df.to_csv(output_dist, sep='\t', index=False)

    print(f"Written threshold evaluation to {output_eval}")
    print(f"Written optimal thresholds to {output_optimal}")
    print(f"Written score distributions to {output_dist}")


if __name__ == '__main__':
    main()
