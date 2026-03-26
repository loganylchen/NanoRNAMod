"""
Generate detailed per-tool benchmark reports with site-by-site comparison.

Output includes:
- For each predicted site: transcript, position, score, truth_match, distance_to_truth
- For each truth site: transcript, position, modification_type, was_detected, nearest_prediction_distance
"""

import os
import sys
import pandas as pd
import numpy as np

# Import shared utilities
# When Snakemake copies script to temp location, we need to find original scripts dir
# Strategy: Check multiple possible locations for benchmark_utils.py

def find_benchmark_utils():
    """Find the directory containing benchmark_utils.py."""
    # 1. Check script's own directory (direct execution or when copied together)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if os.path.exists(os.path.join(script_dir, "benchmark_utils.py")):
        return script_dir

    # 2. Check relative to current working directory
    cwd = os.getcwd()
    candidates = [
        os.path.join(cwd, "workflow", "scripts"),
        os.path.join(cwd, "scripts"),
        cwd,  # If running directly from workflow/scripts
    ]
    for candidate in candidates:
        if os.path.exists(os.path.join(candidate, "benchmark_utils.py")):
            return candidate

    # 3. Search upward from cwd for workflow/scripts structure
    search_dir = cwd
    for _ in range(5):  # Max 5 levels up
        scripts_path = os.path.join(search_dir, "workflow", "scripts")
        if os.path.exists(os.path.join(scripts_path, "benchmark_utils.py")):
            return scripts_path
        parent = os.path.dirname(search_dir)
        if parent == search_dir:
            break
        search_dir = parent

    # 4. Check if it's a sibling of the temp script (same .snakemake/scripts dir)
    if ".snakemake" in script_dir:
        # Try to find original workflow/scripts from snakemake metadata
        # The original Snakefile location might be available
        pass

    return None

_utils_dir = find_benchmark_utils()
if _utils_dir and _utils_dir not in sys.path:
    sys.path.insert(0, _utils_dir)

from benchmark_utils import tool_from_path, normalize_columns


def detect_score_column(df, tool_name=None):
    """Detect the most likely score column for ranking predictions."""
    tool_score_map = {
        'xpore': ['p_value', 'diff_mod', 'diff_mod_frac', 'mod_ratio'],
        'nanocompore': ['pvalue', 'logit_pvalue', 'logit', 'p_value'],
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


def find_nearest_truth(pred_row, truth_df, window=0):
    """
    Find the nearest truth site to a prediction.

    Returns:
        tuple: (found, distance, truth_mod_type)
    """
    tx = pred_row['transcript']
    pos = pred_row['position']

    # Filter to same transcript
    truth_tx = truth_df[truth_df['transcript'] == tx]

    if truth_tx.empty:
        return (False, None, None)

    # Calculate distances
    distances = np.abs(truth_tx['position'].values - pos)
    min_idx = np.argmin(distances)
    min_dist = distances[min_idx]

    if min_dist <= window:
        nearest_row = truth_tx.iloc[min_idx]
        mod_type = nearest_row.get('modification_type', 'unknown')
        return (True, min_dist, mod_type)

    return (False, min_dist, None)


def find_nearest_prediction(truth_row, pred_df, window=0):
    """
    Find the nearest prediction to a truth site.

    Returns:
        tuple: (found, distance, pred_score)
    """
    tx = truth_row['transcript']
    pos = truth_row['position']

    # Filter to same transcript
    pred_tx = pred_df[pred_df['transcript'] == tx]

    if pred_tx.empty:
        return (False, None, None)

    # Calculate distances
    distances = np.abs(pred_tx['position'].values - pos)
    min_idx = np.argmin(distances)
    min_dist = distances[min_idx]

    if min_dist <= window:
        nearest_row = pred_tx.iloc[min_idx]
        score = nearest_row.get('score', None)
        return (True, min_dist, score)

    return (False, min_dist, None)


def generate_detailed_report(pred_df, truth_pos, tool_name, score_col, window=0):
    """
    Generate detailed site-by-site comparison for a tool.

    Returns:
        tuple: (predictions_report, truth_report)
    """
    # Normalize score column name
    if score_col and score_col in pred_df.columns:
        pred_df = pred_df.copy()
        pred_df['score'] = pred_df[score_col]
    elif score_col:
        pred_df = pred_df.copy()
        pred_df['score'] = np.nan
    else:
        pred_df = pred_df.copy()
        pred_df['score'] = np.nan

    # Report for predictions
    pred_records = []
    for _, pred_row in pred_df.iterrows():
        found, dist, mod_type = find_nearest_truth(pred_row, truth_pos, window)

        pred_records.append({
            'tool': tool_name,
            'transcript': pred_row['transcript'],
            'position': int(pred_row['position']),
            'score': pred_row['score'],
            'truth_match': 'yes' if found else 'no',
            'distance_to_truth': dist if found else None,
            'matched_modification': mod_type if found else None,
        })

    pred_report = pd.DataFrame(pred_records)

    # Report for truth sites
    truth_records = []
    for _, truth_row in truth_pos.iterrows():
        found, dist, score = find_nearest_prediction(truth_row, pred_df, window)

        truth_records.append({
            'tool': tool_name,
            'transcript': truth_row['transcript'],
            'position': int(truth_row['position']),
            'modification_type': truth_row.get('modification_type', 'unknown'),
            'detected': 'yes' if found else 'no',
            'nearest_prediction_distance': dist,
            'nearest_prediction_score': score,
        })

    truth_report = pd.DataFrame(truth_records)

    return pred_report, truth_report


def main():
    # When called by Snakemake
    if 'snakemake' in globals():
        result_files = snakemake.input.results
        truth_set_path = snakemake.input.truth_set
        output_pred = snakemake.output.predictions
        output_truth = snakemake.output.truth
        window_param = snakemake.params.get('window', 0)

        if isinstance(window_param, list):
            windows = window_param
        else:
            windows = [int(window_param)]
    else:
        import argparse
        parser = argparse.ArgumentParser(
            description='Generate detailed per-tool benchmark reports'
        )
        parser.add_argument('--results', nargs='+', required=True,
                          help='Tool result TSV files')
        parser.add_argument('--truth', required=True,
                          help='Ground truth TSV file')
        parser.add_argument('--output-pred', required=True,
                          help='Output file for predictions report')
        parser.add_argument('--output-truth', required=True,
                          help='Output file for truth report')
        parser.add_argument('--window', type=int, default=0,
                          help='Positional tolerance window')
        args = parser.parse_args()

        result_files = args.results
        truth_set_path = args.truth
        output_pred = args.output_pred
        output_truth = args.output_truth
        windows = [args.window]

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
            'tool', 'transcript', 'position', 'score', 'truth_match',
            'distance_to_truth', 'matched_modification'
        ]).to_csv(output_pred, sep='\t', index=False)
        pd.DataFrame(columns=[
            'tool', 'transcript', 'position', 'modification_type',
            'detected', 'nearest_prediction_distance', 'nearest_prediction_score'
        ]).to_csv(output_truth, sep='\t', index=False)
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
                if tool not in tool_dfs:
                    tool_dfs[tool] = []
                    if score_col:
                        tool_score_cols[tool] = score_col
                tool_dfs[tool].append(df)
        except Exception as e:
            print(f"Warning: could not read {f}: {e}")

    # Concatenate results per tool
    for tool in tool_dfs:
        tool_dfs[tool] = pd.concat(tool_dfs[tool], ignore_index=True).drop_duplicates(
            subset=['transcript', 'position']
        )

    # Generate reports
    all_pred_reports = []
    all_truth_reports = []

    for tool, pred_df in tool_dfs.items():
        score_col = tool_score_cols.get(tool)

        # Note: Detailed reports use the first window only.
        # For multi-window analysis, see accuracy_summary.tsv
        window = windows[0]

        pred_rep, truth_rep = generate_detailed_report(
            pred_df, truth_pos, tool, score_col, window
        )

        all_pred_reports.append(pred_rep)
        all_truth_reports.append(truth_rep)

    # Combine and write
    if all_pred_reports:
        combined_pred = pd.concat(all_pred_reports, ignore_index=True)
        combined_truth = pd.concat(all_truth_reports, ignore_index=True)
    else:
        combined_pred = pd.DataFrame(columns=[
            'tool', 'transcript', 'position', 'score', 'truth_match',
            'distance_to_truth', 'matched_modification'
        ])
        combined_truth = pd.DataFrame(columns=[
            'tool', 'transcript', 'position', 'modification_type',
            'detected', 'nearest_prediction_distance', 'nearest_prediction_score'
        ])

    combined_pred.to_csv(output_pred, sep='\t', index=False)
    combined_truth.to_csv(output_truth, sep='\t', index=False)

    print(f"Written predictions report to {output_pred}")
    print(f"Written truth report to {output_truth}")


if __name__ == '__main__':
    main()
