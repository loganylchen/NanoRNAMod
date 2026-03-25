#!/usr/bin/env python3
"""
Shared utilities for benchmark scripts.

This module contains common functions used across multiple benchmark scripts
to avoid code duplication and ensure consistency.

Author: NanoRNAMod Development Team
"""

import os
import re
import pandas as pd
import numpy as np


def tool_from_path(path):
    """
    Extract tool name from file path.

    Args:
        path: File path containing tool name (e.g., results/modifications/xpore/...)

    Returns:
        Tool name string (e.g., 'xpore')
    """
    parts = path.split(os.sep)
    for i, part in enumerate(parts):
        if part == "modifications" and i + 1 < len(parts):
            return parts[i + 1]
    # Fallback: use parent directory of the file
    return os.path.basename(os.path.dirname(os.path.dirname(path)))


def comparison_from_path(path):
    """
    Extract comparison name from file path.

    Args:
        path: File path containing comparison name

    Returns:
        Comparison name string (e.g., 'A_F') or 'unknown' if not found
    """
    parts = path.split(os.sep)
    for i, part in enumerate(parts):
        if part == "modifications" and i + 2 < len(parts):
            return parts[i + 2]
    return "unknown"


def normalize_columns(df):
    """
    Normalize various column names to standard 'transcript' and 'position'.

    Handles different naming conventions across tools:
    - transcript: id, ref_id, chrom
    - position: pos, start, start_loc, transcript_pos, transcript_loc

    Args:
        df: DataFrame with potentially non-standard column names

    Returns:
        DataFrame with normalized column names
    """
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


def match_column(df_columns, pattern):
    """
    Match a pattern against DataFrame columns using multiple strategies.

    Args:
        df_columns: List of column names
        pattern: Pattern to match (case-insensitive)

    Returns:
        Matched column name or None
    """
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


def is_pvalue_column(col_name):
    """
    Determine if a column is a p-value (lower=better) or score (higher=better).

    Args:
        col_name: Column name to check

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
    if any(x in col_lower for x in ['p_value', 'pvalue', 'pval', 'padj',
                                      'adjusted_p_value', 'fdr', 'q_value', 'qval']):
        return True

    # Probability columns (higher = more modification, better)
    if any(x in col_lower for x in ['probability', 'prob', 'psi_probability', 'mean_p_mod']):
        return False

    # Stoichiometry (higher = more modification)
    if 'stoichiometry' in col_lower:
        return False

    # Effect size (higher = larger difference)
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


def load_truth_set(truth_path):
    """
    Load and validate truth set file.

    Args:
        truth_path: Path to truth set TSV file

    Returns:
        Tuple of (truth_df, truth_pos_df)
    """
    truth_df = pd.read_csv(truth_path, sep='\t')

    if 'label' in truth_df.columns:
        truth_pos = truth_df[truth_df['label'] != '-'].copy()
    else:
        truth_pos = truth_df.copy()

    truth_pos = normalize_columns(truth_pos)

    return truth_df, truth_pos


def compute_f1(precision, recall):
    """
    Compute F1 score from precision and recall.

    Args:
        precision: Precision value
        recall: Recall value

    Returns:
        F1 score (harmonic mean of precision and recall)
    """
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


def bin_coverage(coverage_values, bins=None):
    """
    Bin coverage values into discrete categories.

    Args:
        coverage_values: Array of coverage depths
        bins: List of bin boundaries (default: [0, 10, 20, 50, 100, 200, 500])

    Returns:
        Array of bin labels as strings
    """
    if bins is None:
        bins = [0, 10, 20, 50, 100, 200, 500]

    bin_edges = bins + [np.inf]

    labels = []
    for i in range(len(bin_edges) - 1):
        if bin_edges[i + 1] == np.inf:
            labels.append(f"[{bin_edges[i]},inf)")
        else:
            labels.append(f"[{bin_edges[i]},{bin_edges[i+1]})")

    indices = np.digitize(coverage_values, bin_edges[1:-1], right=False)
    result = np.array([labels[min(i, len(labels)-1)] for i in indices])

    return result


def generate_thresholds(scores, n_thresholds=100):
    """
    Generate threshold values based on score distribution.

    Uses percentile-based thresholds for even coverage across the score range.

    Args:
        scores: Series of score values
        n_thresholds: Number of thresholds to generate

    Returns:
        Sorted list of unique threshold values
    """
    scores = scores.dropna()

    if scores.empty:
        return []

    percentiles = np.linspace(0, 100, n_thresholds)
    percentile_thresholds = np.percentile(scores, percentiles)

    thresholds = sorted(set(percentile_thresholds))

    return thresholds
