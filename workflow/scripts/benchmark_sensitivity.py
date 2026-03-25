#!/usr/bin/env python3
"""
Benchmark Sensitivity Analysis Module

Analyzes how benchmark metrics vary with:
1. Coverage depth (accuracy vs read depth)
2. Score distribution (separation between TP and FP)
3. Threshold robustness (stability across data splits)

Outputs:
    - coverage_analysis.tsv: Metrics by coverage depth bin
    - score_distribution.tsv: Score distributions by truth status
    - threshold_robustness.tsv: Threshold stability across CV splits

Author: NanoRNAMod Development Team
Date: 2026-03-25
"""

import os
import sys
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import ks_2samp, mannwhitneyu
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from tqdm import tqdm
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# Coverage Depth Analysis
# =============================================================================

def bin_coverage(
    coverage_values: np.ndarray,
    bins: List[int] = None
) -> np.ndarray:
    """
    Bin coverage values into discrete categories.

    Args:
        coverage_values: Array of coverage depths
        bins: List of bin boundaries (e.g., [0, 10, 20, 50, 100, 200, 500])

    Returns:
        Array of bin labels as strings
    """
    if bins is None:
        bins = [0, 10, 20, 50, 100, 200, 500]

    # Add infinity as upper bound
    bin_edges = bins + [np.inf]

    labels = []
    for i in range(len(bin_edges) - 1):
        if bin_edges[i + 1] == np.inf:
            labels.append(f"[{bin_edges[i]},inf)")
        else:
            labels.append(f"[{bin_edges[i]},{bin_edges[i+1]})")

    # Use numpy digitize for binning
    indices = np.digitize(coverage_values, bin_edges[1:-1], right=False)
    result = np.array([labels[min(i, len(labels)-1)] for i in indices])

    return result


def compute_metrics_at_threshold(
    y_true: np.ndarray,
    y_scores: np.ndarray,
    threshold: float
) -> Dict[str, float]:
    """
    Compute precision, recall, F1 at a specific threshold.

    Args:
        y_true: Binary truth labels
        y_scores: Prediction scores
        threshold: Classification threshold

    Returns:
        Dictionary with precision, recall, f1, specificity
    """
    y_pred = (y_scores >= threshold).astype(int)

    TP = np.sum((y_pred == 1) & (y_true == 1))
    FP = np.sum((y_pred == 1) & (y_true == 0))
    TN = np.sum((y_pred == 0) & (y_true == 0))
    FN = np.sum((y_pred == 0) & (y_true == 1))

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0.0

    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'specificity': specificity
    }


def analyze_coverage_effect(
    predictions_df: pd.DataFrame,
    coverage_bins: List[int] = None,
    threshold: float = 0.5
) -> pd.DataFrame:
    """
    Analyze how metrics vary with coverage depth.

    Args:
        predictions_df: DataFrame with columns [tool, coverage, score, truth, comparison]
        coverage_bins: Coverage depth bin boundaries
        threshold: Score threshold for classification

    Returns:
        DataFrame with metrics by coverage bin and tool
    """
    if coverage_bins is None:
        coverage_bins = [0, 10, 20, 50, 100, 200, 500]

    # Ensure coverage column exists
    if 'coverage' not in predictions_df.columns:
        logger.warning("No coverage column found, using 'n_reads' as proxy")
        if 'n_reads' in predictions_df.columns:
            predictions_df['coverage'] = predictions_df['n_reads']
        else:
            logger.error("No coverage information available")
            return pd.DataFrame()

    results = []
    tools = predictions_df['tool'].unique()

    for tool in tools:
        tool_df = predictions_df[predictions_df['tool'] == tool].copy()

        # Bin coverage
        tool_df['coverage_bin'] = bin_coverage(
            tool_df['coverage'].values,
            bins=coverage_bins
        )

        for bin_label in sorted(tool_df['coverage_bin'].unique()):
            bin_df = tool_df[tool_df['coverage_bin'] == bin_label]

            if len(bin_df) < 5:
                continue

            y_true = bin_df['truth'].values
            y_scores = bin_df['score'].values

            metrics = compute_metrics_at_threshold(y_true, y_scores, threshold)

            results.append({
                'tool': tool,
                'coverage_bin': bin_label,
                'n_sites': len(bin_df),
                'n_positive': int(y_true.sum()),
                'n_negative': int((1 - y_true).sum()),
                'mean_coverage': bin_df['coverage'].mean(),
                **metrics
            })

    return pd.DataFrame(results)


# =============================================================================
# Score Distribution Analysis
# =============================================================================

def analyze_score_distribution(
    predictions_df: pd.DataFrame,
    tools: List[str] = None
) -> pd.DataFrame:
    """
    Analyze score distributions by truth status.

    Computes distribution statistics and KS tests to quantify
    separation between true positive and false positive scores.

    Args:
        predictions_df: DataFrame with columns [tool, score, truth]
        tools: List of tools to analyze (default: all)

    Returns:
        DataFrame with distribution statistics and test results
    """
    if tools is None:
        tools = predictions_df['tool'].unique()

    results = []

    for tool in tools:
        tool_df = predictions_df[predictions_df['tool'] == tool]

        if 'score' not in tool_df.columns:
            logger.warning(f"No score column for tool {tool}")
            continue

        # Separate by truth status
        pos_scores = tool_df[tool_df['truth'] == 1]['score'].dropna().values
        neg_scores = tool_df[tool_df['truth'] == 0]['score'].dropna().values

        if len(pos_scores) < 2 or len(neg_scores) < 2:
            continue

        # Distribution statistics for positives
        pos_stats = {
            'pos_mean': np.mean(pos_scores),
            'pos_median': np.median(pos_scores),
            'pos_std': np.std(pos_scores),
            'pos_q25': np.percentile(pos_scores, 25),
            'pos_q75': np.percentile(pos_scores, 75),
            'pos_min': np.min(pos_scores),
            'pos_max': np.max(pos_scores),
            'n_pos': len(pos_scores)
        }

        # Distribution statistics for negatives
        neg_stats = {
            'neg_mean': np.mean(neg_scores),
            'neg_median': np.median(neg_scores),
            'neg_std': np.std(neg_scores),
            'neg_q25': np.percentile(neg_scores, 25),
            'neg_q75': np.percentile(neg_scores, 75),
            'neg_min': np.min(neg_scores),
            'neg_max': np.max(neg_scores),
            'n_neg': len(neg_scores)
        }

        # Separation metrics
        mean_diff = pos_stats['pos_mean'] - neg_stats['neg_mean']
        pooled_std = np.sqrt(
            (pos_stats['pos_std']**2 + neg_stats['neg_std']**2) / 2
        )
        cohens_d = mean_diff / pooled_std if pooled_std > 0 else 0.0

        # Overlap coefficient (histogram-based)
        n_bins = 50
        pos_hist, bin_edges = np.histogram(pos_scores, bins=n_bins, density=True)
        neg_hist, _ = np.histogram(neg_scores, bins=bin_edges, density=True)
        overlap = np.minimum(pos_hist, neg_hist).sum() / max(pos_hist.sum(), neg_hist.sum())

        # KS test
        ks_stat, ks_pvalue = ks_2samp(pos_scores, neg_scores)

        # Mann-Whitney U test
        mw_stat, mw_pvalue = mannwhitneyu(pos_scores, neg_scores, alternative='two-sided')

        results.append({
            'tool': tool,
            'mean_diff': mean_diff,
            'cohens_d': cohens_d,
            'overlap_coefficient': overlap,
            'ks_statistic': ks_stat,
            'ks_pvalue': ks_pvalue,
            'mannwhitney_stat': mw_stat,
            'mannwhitney_pvalue': mw_pvalue,
            **pos_stats,
            **neg_stats
        })

    return pd.DataFrame(results)


# =============================================================================
# Threshold Robustness Analysis
# =============================================================================

def find_optimal_threshold(
    y_true: np.ndarray,
    y_scores: np.ndarray,
    criterion: str = 'f1'
) -> Tuple[float, float]:
    """
    Find optimal threshold based on criterion.

    Args:
        y_true: Binary truth labels
        y_scores: Prediction scores
        criterion: Optimization criterion ('f1', 'youden', 'precision')

    Returns:
        Tuple of (optimal_threshold, best_metric_value)
    """
    if len(np.unique(y_true)) < 2:
        return 0.5, 0.0

    if criterion == 'f1':
        precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
        f1_scores = 2 * precision * recall / (precision + recall + 1e-10)
        best_idx = np.argmax(f1_scores)
        return thresholds[best_idx] if best_idx < len(thresholds) else thresholds[-1], f1_scores[best_idx]

    elif criterion == 'youden':
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        youden = tpr - fpr
        best_idx = np.argmax(youden)
        return thresholds[best_idx], youden[best_idx]

    elif criterion == 'precision':
        precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
        best_idx = np.argmax(precision)
        return thresholds[best_idx] if best_idx < len(thresholds) else thresholds[-1], precision[best_idx]

    return 0.5, 0.0


def analyze_threshold_robustness(
    predictions_df: pd.DataFrame,
    n_splits: int = 5,
    criterion: str = 'f1',
    seed: int = 42
) -> pd.DataFrame:
    """
    Assess threshold stability via cross-validation.

    Args:
        predictions_df: DataFrame with columns [tool, score, truth, comparison]
        n_splits: Number of CV splits
        criterion: Optimization criterion
        seed: Random seed

    Returns:
        DataFrame with threshold variance analysis
    """
    np.random.seed(seed)

    results = []
    tools = predictions_df['tool'].unique()

    for tool in tools:
        tool_df = predictions_df[predictions_df['tool'] == tool].copy()

        if 'score' not in tool_df.columns:
            continue

        # Remove NaN scores
        tool_df = tool_df.dropna(subset=['score'])

        if len(tool_df) < n_splits * 10:
            continue

        y_true = tool_df['truth'].values
        y_scores = tool_df['score'].values

        # Cross-validation
        thresholds = []
        metrics = []

        try:
            skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)

            for train_idx, _ in skf.split(y_scores, y_true):
                train_true = y_true[train_idx]
                train_scores = y_scores[train_idx]

                thresh, metric = find_optimal_threshold(
                    train_true, train_scores, criterion=criterion
                )
                thresholds.append(thresh)
                metrics.append(metric)

            thresholds = np.array(thresholds)
            metrics = np.array(metrics)

            # Compute stability metrics
            threshold_mean = np.mean(thresholds)
            threshold_std = np.std(thresholds)
            threshold_cv = threshold_std / threshold_mean if threshold_mean > 0 else 0.0
            threshold_range = np.max(thresholds) - np.min(thresholds)

            metric_mean = np.mean(metrics)
            metric_std = np.std(metrics)

            # Full-data optimal threshold
            full_thresh, full_metric = find_optimal_threshold(
                y_true, y_scores, criterion=criterion
            )

            results.append({
                'tool': tool,
                'n_splits': n_splits,
                'criterion': criterion,
                'threshold_mean': threshold_mean,
                'threshold_std': threshold_std,
                'threshold_cv': threshold_cv,
                'threshold_min': np.min(thresholds),
                'threshold_max': np.max(thresholds),
                'threshold_range': threshold_range,
                'metric_mean': metric_mean,
                'metric_std': metric_std,
                'full_data_threshold': full_thresh,
                'full_data_metric': full_metric,
                'n_samples': len(tool_df)
            })

        except Exception as e:
            logger.warning(f"Could not compute threshold robustness for {tool}: {e}")
            continue

    return pd.DataFrame(results)


# =============================================================================
# Main Entry Point
# =============================================================================

def main():
    """Main entry point for sensitivity analysis."""

    # Parse Snakemake parameters if available
    try:
        predictions_input = snakemake.input.predictions
        truth_set_input = snakemake.input.truth_set

        output_coverage = snakemake.output.coverage
        output_score_dist = snakemake.output.score_dist
        output_threshold_robust = snakemake.output.threshold_robust

        coverage_bins = snakemake.params.get('coverage_bins', [0, 10, 20, 50, 100, 200, 500])
        n_splits = snakemake.params.get('n_splits', 5)

        logger.info("Running in Snakemake mode")

    except NameError:
        import argparse

        parser = argparse.ArgumentParser(description='Benchmark Sensitivity Analysis')
        parser.add_argument('--predictions', required=True, help='Detailed predictions TSV')
        parser.add_argument('--truth-set', required=True, help='Truth set TSV')
        parser.add_argument('--output-coverage', required=True, help='Coverage analysis output')
        parser.add_argument('--output-score-dist', required=True, help='Score distribution output')
        parser.add_argument('--output-threshold-robust', required=True, help='Threshold robustness output')
        parser.add_argument('--coverage-bins', type=int, nargs='+',
                            default=[0, 10, 20, 50, 100, 200, 500])
        parser.add_argument('--n-splits', type=int, default=5)

        args = parser.parse_args()

        predictions_input = args.predictions
        truth_set_input = args.truth_set
        output_coverage = args.output_coverage
        output_score_dist = args.output_score_dist
        output_threshold_robust = args.output_threshold_robust
        coverage_bins = args.coverage_bins
        n_splits = args.n_splits

        logger.info("Running in command line mode")

    # Create output directories
    for output_path in [output_coverage, output_score_dist, output_threshold_robust]:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Load data
    logger.info(f"Loading predictions from: {predictions_input}")
    predictions_df = pd.read_csv(predictions_input, sep='\t')

    logger.info(f"Loading truth set from: {truth_set_input}")
    truth_df = pd.read_csv(truth_set_input, sep='\t')

    # Ensure required columns
    required_cols = ['tool', 'score', 'truth']
    for col in required_cols:
        if col not in predictions_df.columns:
            logger.error(f"Required column '{col}' not found in predictions")
            sys.exit(1)

    # =========================================================================
    # Coverage Analysis
    # =========================================================================
    logger.info("=" * 60)
    logger.info("PHASE 1: Coverage Depth Analysis")
    logger.info("=" * 60)

    coverage_df = analyze_coverage_effect(
        predictions_df,
        coverage_bins=coverage_bins,
        threshold=0.5
    )

    if not coverage_df.empty:
        coverage_df.to_csv(output_coverage, sep='\t', index=False)
        logger.info(f"Saved coverage analysis to: {output_coverage}")
        logger.info(f"  - {len(coverage_df)} bin-tool combinations analyzed")
    else:
        logger.warning("Coverage analysis produced no results")
        pd.DataFrame().to_csv(output_coverage, sep='\t', index=False)

    # =========================================================================
    # Score Distribution Analysis
    # =========================================================================
    logger.info("=" * 60)
    logger.info("PHASE 2: Score Distribution Analysis")
    logger.info("=" * 60)

    score_dist_df = analyze_score_distribution(predictions_df)

    if not score_dist_df.empty:
        score_dist_df.to_csv(output_score_dist, sep='\t', index=False)
        logger.info(f"Saved score distributions to: {output_score_dist}")
        logger.info(f"  - {len(score_dist_df)} tools analyzed")

        # Log key findings
        for _, row in score_dist_df.iterrows():
            logger.info(
                f"  {row['tool']}: Cohen's d={row['cohens_d']:.2f}, "
                f"KS={row['ks_statistic']:.3f}, Overlap={row['overlap_coefficient']:.2f}"
            )
    else:
        logger.warning("Score distribution analysis produced no results")
        pd.DataFrame().to_csv(output_score_dist, sep='\t', index=False)

    # =========================================================================
    # Threshold Robustness Analysis
    # =========================================================================
    logger.info("=" * 60)
    logger.info("PHASE 3: Threshold Robustness Analysis")
    logger.info("=" * 60)

    threshold_df = analyze_threshold_robustness(
        predictions_df,
        n_splits=n_splits,
        criterion='f1'
    )

    if not threshold_df.empty:
        threshold_df.to_csv(output_threshold_robust, sep='\t', index=False)
        logger.info(f"Saved threshold robustness to: {output_threshold_robust}")
        logger.info(f"  - {len(threshold_df)} tools analyzed")

        # Log key findings
        for _, row in threshold_df.iterrows():
            logger.info(
                f"  {row['tool']}: threshold={row['threshold_mean']:.3f} ± {row['threshold_std']:.3f} "
                f"(CV={row['threshold_cv']:.2%})"
            )
    else:
        logger.warning("Threshold robustness analysis produced no results")
        pd.DataFrame().to_csv(output_threshold_robust, sep='\t', index=False)

    logger.info("=" * 60)
    logger.info("Sensitivity analysis complete!")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
