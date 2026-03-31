"""
Modification ratio extraction for benchmarking.

Extracts modification ratio (0-1) from tools that report it, merges with
ground truth, and outputs a TSV with native and fair evaluation modes.

Tools and their mod_ratio columns:
  - pybaleen: mod_ratio
  - xpore: diff_mod_rate_CASE_vs_CONTROL
  - tandemmod: stoichiometry
  - m6atm: stoichiometry
  - directrm: probability
  - rnano: probability
  - nanopsu: psi_probability
  - nanomud: probability
  - penguin: psi_probability
  - eligos2: error_rate
  - nanocompore: parsed from cluster_counts (GMM_n_clust == 2)
  - baleen: native_mod_ratio (if present)
"""

import sys
import re
import logging
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
)
logger = logging.getLogger(__name__)

# Tool -> column name that represents modification ratio (0-1)
TOOL_MOD_RATIO_COL = {
    'pybaleen': 'mod_ratio',
    'xpore': 'diff_mod_rate_CASE_vs_CONTROL',
    'tandemmod': 'stoichiometry',
    'm6atm': 'stoichiometry',
    'directrm': 'probability',
    'rnano': 'probability',
    'nanopsu': 'psi_probability',
    'nanomud': 'probability',
    'penguin': 'psi_probability',
    'eligos2': 'error_rate',
    'nanocompore': '__parse_cluster_counts__',
    'baleen': 'native_mod_ratio',
}

# Column normalization: various tool column names -> transcript, position
TRANSCRIPT_COLS = ['transcript', 'id', 'ref_id', 'chrom', 'transcript_id', 'contig']
POSITION_COLS = ['position', 'pos', 'start', 'start_loc', 'transcript_pos', 'transcript_loc']


def normalize_columns(df):
    """Normalize transcript and position column names."""
    df = df.copy()
    for col in TRANSCRIPT_COLS:
        if col in df.columns and col != 'transcript':
            df = df.rename(columns={col: 'transcript'})
            break
    for col in POSITION_COLS:
        if col in df.columns and col != 'position':
            df = df.rename(columns={col: 'position'})
            break
    return df


def parse_nanocompore_mod_ratio(df):
    """Parse modification ratio from nanocompore cluster_counts.

    When GMM_n_clust == 2, cluster_counts looks like:
        Control_1:14/18__Native_1:8/22
    The native mod ratio is the fraction of reads in the minor cluster
    (assumed modified) for the native sample.

    When GMM_n_clust == 1, mod_ratio = 0 (no evidence of modification).
    """
    if 'cluster_counts' not in df.columns or 'GMM_n_clust' not in df.columns:
        logger.warning("nanocompore: missing cluster_counts or GMM_n_clust columns")
        return pd.Series(np.nan, index=df.index)

    ratios = pd.Series(0.0, index=df.index)

    mask_2 = df['GMM_n_clust'] == 2
    for idx in df.index[mask_2]:
        cc = str(df.loc[idx, 'cluster_counts'])
        try:
            parts = cc.split('__')
            native_parts = [p for p in parts if 'Native' in p or 'native' in p.lower()]
            if not native_parts:
                ratios[idx] = np.nan
                continue

            # Each native part: Native_1:8/22
            native_fracs = []
            for np_ in native_parts:
                match = re.search(r':(\d+)/(\d+)', np_)
                if match:
                    num, denom = int(match.group(1)), int(match.group(2))
                    if denom > 0:
                        native_fracs.append(num / denom)

            if len(native_fracs) == 0:
                ratios[idx] = np.nan
            elif len(native_fracs) == 1:
                frac = native_fracs[0]
                ratios[idx] = min(frac, 1 - frac)
            else:
                ratios[idx] = min(native_fracs)
        except Exception:
            ratios[idx] = np.nan

    return ratios


def tool_from_path(path):
    """Extract tool name from result file path."""
    import os
    parts = os.path.normpath(path).split(os.sep)
    for i, p in enumerate(parts):
        if p == 'modifications' and i + 1 < len(parts):
            return parts[i + 1]
    basename = os.path.basename(path)
    return basename.replace('_results.tsv', '')


def extract_mod_ratio(df, tool_name):
    """Extract modification ratio column for a tool. Returns Series or None."""
    if tool_name not in TOOL_MOD_RATIO_COL:
        return None

    col = TOOL_MOD_RATIO_COL[tool_name]

    if col == '__parse_cluster_counts__':
        return parse_nanocompore_mod_ratio(df)

    if col not in df.columns:
        logger.warning(f"{tool_name}: column '{col}' not found in results")
        return None

    return pd.to_numeric(df[col], errors='coerce')


def load_truth_set(truth_path, window=0):
    """Load truth set and return set of (transcript, position) tuples."""
    truth = pd.read_csv(truth_path, sep='\t')
    truth = normalize_columns(truth)
    if 'label' in truth.columns:
        truth = truth[truth['label'] != '-']
    positions = set()
    for _, row in truth.iterrows():
        t = row['transcript']
        p = int(row['position'])
        for w in range(-window, window + 1):
            positions.add((t, p + w))
    return positions


def build_mod_ratio_data(result_paths, truth_path, union_predictions_path, window=0):
    """Build modification ratio dataframe for all tools.

    Returns DataFrame with columns:
        tool, comparison, transcript, position, mod_ratio, truth, mode
    where mode is 'native' or 'fair'.
    """
    import os
    truth_positions = load_truth_set(truth_path, window)

    # Load union predictions for fair mode
    union_sites = {}
    if union_predictions_path:
        try:
            union_df = pd.read_csv(union_predictions_path, sep='\t')
            union_df = normalize_columns(union_df)
            for _, row in union_df.iterrows():
                key = (str(row['transcript']), int(row['position']))
                union_sites[key] = True
        except Exception as e:
            logger.warning(f"Could not load union predictions: {e}")

    all_records = []

    for path in result_paths:
        tool = tool_from_path(path)
        if tool not in TOOL_MOD_RATIO_COL:
            continue

        try:
            df = pd.read_csv(path, sep='\t')
        except Exception as e:
            logger.warning(f"Could not read {path}: {e}")
            continue

        df = normalize_columns(df)
        if 'transcript' not in df.columns or 'position' not in df.columns:
            logger.warning(f"{tool}: missing transcript/position columns")
            continue

        mod_ratios = extract_mod_ratio(df, tool)
        if mod_ratios is None:
            continue

        df = df.copy()
        df['mod_ratio'] = mod_ratios
        df['position'] = df['position'].astype(int)

        # Extract comparison from path
        parts = os.path.normpath(path).split(os.sep)
        comparison = 'unknown'
        for i, p in enumerate(parts):
            if p == 'modifications' and i + 2 < len(parts):
                comparison = parts[i + 2]
                break

        # ---- Native mode ----
        # Only sites this tool actually called
        tool_mr_lookup = {}
        for _, row in df.iterrows():
            key = (str(row['transcript']), int(row['position']))
            is_positive = key in truth_positions
            mr = row['mod_ratio']
            if pd.notna(mr):
                tool_mr_lookup[key] = float(mr)
                all_records.append({
                    'tool': tool,
                    'comparison': comparison,
                    'transcript': row['transcript'],
                    'position': int(row['position']),
                    'mod_ratio': float(mr),
                    'truth': 'positive' if is_positive else 'negative',
                    'mode': 'native',
                })

        # ---- Fair mode ----
        # All union sites; sites not called by this tool get mod_ratio = 0
        if union_sites:
            for (t, p) in union_sites:
                is_positive = (t, p) in truth_positions
                mr = tool_mr_lookup.get((t, p), 0.0)
                all_records.append({
                    'tool': tool,
                    'comparison': comparison,
                    'transcript': t,
                    'position': p,
                    'mod_ratio': mr,
                    'truth': 'positive' if is_positive else 'negative',
                    'mode': 'fair',
                })

    return pd.DataFrame(all_records)


def main():
    if 'snakemake' in globals():
        result_files = snakemake.input.results
        truth_set_path = snakemake.input.truth_set
        union_predictions_path = snakemake.input.get('union_predictions', None)
        output_data = snakemake.output.data
        window_param = snakemake.params.get('window', 0)
        if isinstance(window_param, list):
            window = int(window_param[0])
        else:
            window = int(window_param)
        log_file = snakemake.log[0]

        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
        logger.setLevel(logging.INFO)
    else:
        import argparse
        parser = argparse.ArgumentParser(description='Modification ratio extraction')
        parser.add_argument('--results', nargs='+', required=True)
        parser.add_argument('--truth', required=True)
        parser.add_argument('--union-predictions', default=None)
        parser.add_argument('--output-data', required=True)
        parser.add_argument('--window', type=int, default=0)
        args = parser.parse_args()

        result_files = args.results
        truth_set_path = args.truth
        union_predictions_path = args.union_predictions
        output_data = args.output_data
        window = args.window

    if isinstance(result_files, str):
        result_files = [result_files]

    logger.info(f"Processing {len(result_files)} result files")
    logger.info(f"Truth set: {truth_set_path}")

    mod_ratio_df = build_mod_ratio_data(
        result_files, truth_set_path, union_predictions_path, window
    )

    if mod_ratio_df.empty:
        logger.warning("No modification ratio data found")
        pd.DataFrame(columns=[
            'tool', 'comparison', 'transcript', 'position',
            'mod_ratio', 'truth', 'mode'
        ]).to_csv(output_data, sep='\t', index=False)
        return

    logger.info(f"Total records: {len(mod_ratio_df)}")
    for tool in mod_ratio_df['tool'].unique():
        t_df = mod_ratio_df[mod_ratio_df['tool'] == tool]
        n_native = len(t_df[t_df['mode'] == 'native'])
        n_fair = len(t_df[t_df['mode'] == 'fair'])
        logger.info(f"  {tool}: {n_native} native, {n_fair} fair")

    mod_ratio_df.to_csv(output_data, sep='\t', index=False)
    logger.info(f"Saved data to {output_data}")


if __name__ == '__main__':
    main()
