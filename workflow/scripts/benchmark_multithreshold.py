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
    """Determine if a column is a p-value (lower=better) or score (higher=better)."""
    if not col_name:
        return False

    col_lower = col_name.lower().replace(' ', '_')

    # Explicit p-value indicators
    if any(x in col_lower for x in ['p_value', 'pvalue', 'pval', 'fdr', 'qval']):
        return True

    # Tool-specific heuristics
    if tool_name == 'xpore' and 'p_value' in col_lower:
        return True
    if tool_name == 'nanocompore' and 'pvalue' in col_lower:
        return True
    if tool_name == 'differr' and '-log10' in col_lower:
        return False  # -log10 is higher=better
    if tool_name in ['tandemmod', 'directrm', 'm6atm', 'rnano'] and 'probability' in col_lower:
        return False

    # Default: assume higher is better for score/probability columns
    return False


def generate_thresholds(scores, n_thresholds=50, min_thresholds=10, custom_thresholds=None):
    """
    Generate threshold values based on score distribution.

    Uses percentiles to ensure even coverage across the score range.
    """
    scores = scores.dropna()

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


def compute_score_distribution(pred_df, score_col):
    """Compute summary statistics for score distribution."""
    scores = pred_df[score_col].dropna()

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
        median_score = pred_df[score_col].median()
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
            'optimal_threshold', 'optimal_f1', 'precision', 'recall', 'tp', 'fp', 'fn', 'criterion'
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

    for f in result_files:
        tool = tool_from_path(f)
        try:
            df = pd.read_csv(f, sep='\t')
            df = normalize_columns(df)

            if 'transcript' in df.columns and 'position' in df.columns:
                score_col = detect_score_column(df, tool)
                if score_col:
                    if tool not in tool_dfs:
                        tool_dfs[tool] = []
                        tool_score_cols[tool] = score_col
                        tool_is_pvalue[tool] = is_pvalue_column(score_col, tool)
                    tool_dfs[tool].append(df)
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

            # Generate thresholds from score distribution
            scores = pred_subset[score_col]
            thresholds = generate_thresholds(scores, n_thresholds=n_thresholds,
                                            custom_thresholds=custom_thresholds)

            for window in windows:
                eval_records, optimal_record, dist_record = evaluate_tool_at_thresholds(
                    pred_subset, truth_subset, tool, mod_type,
                    score_col, is_pvalue, thresholds, window
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
            'threshold', 'tp', 'fp', 'fn', 'precision', 'recall', 'f1', 'predictions'
        ])

    if all_optimal_records:
        optimal_df = pd.DataFrame(all_optimal_records)
        optimal_df = optimal_df.sort_values(
            ['modification_type', 'window', 'optimal_f1'],
            ascending=[True, True, False]
        )
    else:
        optimal_df = pd.DataFrame(columns=[
            'tool', 'modification_type', 'window', 'score_column', 'is_pvalue',
            'optimal_threshold', 'optimal_f1', 'precision', 'recall', 'tp', 'fp', 'fn', 'criterion'
        ])

    if all_dist_records:
        dist_df = pd.DataFrame(all_dist_records)
        dist_df = dist_df.sort_values(['tool', 'modification_type'])
    else:
        dist_df = pd.DataFrame(columns=[
            'tool', 'modification_type', 'score_column', 'is_pvalue',
            'count', 'min', 'max', 'mean', 'std', 'median', 'q25', 'q75'
        ])

    eval_df.to_csv(output_eval, sep='\t', index=False)
    optimal_df.to_csv(output_optimal, sep='\t', index=False)
    dist_df.to_csv(output_dist, sep='\t', index=False)

    print(f"Written threshold evaluation to {output_eval}")
    print(f"Written optimal thresholds to {output_optimal}")
    print(f"Written score distributions to {output_dist}")


if __name__ == '__main__':
    main()
