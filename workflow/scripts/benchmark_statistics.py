#!/usr/bin/env python3
"""
Benchmark Statistical Analysis Module

Computes bootstrap confidence intervals, significance tests, and effect sizes
for RNA modification detection tool benchmarks.

Outputs:
    - bootstrap_ci.tsv: 95% CIs for all metrics by tool and modification type
    - significance_tests.tsv: Paired Wilcoxon and permutation tests
    - fdr_corrected.tsv: FDR-adjusted p-values
    - effect_sizes.tsv: Cohen's d and odds ratios

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
from scipy.stats import wilcoxon, mannwhitneyu, ks_2samp
from tqdm import tqdm
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=RuntimeWarning)


# =============================================================================
# Native FDR Correction (replaces statsmodels dependency)
# =============================================================================

def benjamini_hochberg_correction(p_values: np.ndarray, alpha: float = 0.05) -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply Benjamini-Hochberg FDR correction natively.

    Args:
        p_values: Array of p-values (may contain NaN)
        alpha: Significance level

    Returns:
        Tuple of (rejected, adjusted_p_values)
    """
    p_values = np.asarray(p_values)
    n = len(p_values)

    if n == 0:
        return np.array([], dtype=bool), np.array([])

    # Handle NaN values
    valid_mask = ~np.isnan(p_values)
    valid_p = p_values[valid_mask].copy()

    if len(valid_p) == 0:
        return np.zeros(n, dtype=bool), p_values.copy()

    # Sort p-values and get ranks
    sorted_indices = np.argsort(valid_p)
    sorted_p = valid_p[sorted_indices]
    ranks = np.arange(1, len(sorted_p) + 1)

    # BH formula: adjusted p = p * n / rank
    adjusted = sorted_p * n / ranks

    # Ensure monotonicity (cumulative minimum from right to left)
    for i in range(len(adjusted) - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])

    # Cap at 1.0
    adjusted = np.minimum(adjusted, 1.0)

    # Map back to original order
    full_adj = np.full(n, np.nan)
    full_adj[valid_mask] = adjusted[np.argsort(sorted_indices)]

    # Determine significance
    rejected = full_adj <= alpha

    return rejected, full_adj


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# Bootstrap Confidence Intervals
# =============================================================================

def compute_bootstrap_ci(
    values: np.ndarray,
    n_bootstrap: int = 1000,
    ci: float = 95,
    seed: Optional[int] = None
) -> Tuple[float, float, float]:
    """
    Compute bootstrap confidence interval for a set of values.

    Args:
        values: Array of metric values
        n_bootstrap: Number of bootstrap iterations
        ci: Confidence interval width (default: 95)
        seed: Random seed for reproducibility

    Returns:
        Tuple of (lower_bound, upper_bound, point_estimate)
    """
    if len(values) == 0:
        return np.nan, np.nan, np.nan

    if seed is not None:
        np.random.seed(seed)

    # Point estimate (mean)
    point_estimate = np.nanmean(values)

    # Handle edge cases
    if np.all(np.isnan(values)):
        return np.nan, np.nan, np.nan

    valid_values = values[~np.isnan(values)]
    if len(valid_values) < 2:
        return point_estimate, point_estimate, point_estimate

    # Bootstrap resampling
    bootstrap_means = []
    for _ in range(n_bootstrap):
        sample = np.random.choice(valid_values, size=len(valid_values), replace=True)
        bootstrap_means.append(np.nanmean(sample))

    bootstrap_means = np.array(bootstrap_means)

    # Compute CI percentiles
    alpha = (100 - ci) / 2
    lower = np.percentile(bootstrap_means, alpha)
    upper = np.percentile(bootstrap_means, 100 - alpha)

    return lower, upper, point_estimate


def compute_bootstrap_ci_grouped(
    df: pd.DataFrame,
    metric_col: str,
    group_cols: List[str],
    n_bootstrap: int = 1000,
    ci: float = 95,
    seed: int = 42
) -> pd.DataFrame:
    """
    Compute bootstrap CIs for a metric stratified by groups.

    Args:
        df: DataFrame with metric values
        metric_col: Column name for the metric
        group_cols: Columns to stratify by
        n_bootstrap: Number of bootstrap iterations
        ci: Confidence interval width
        seed: Random seed

    Returns:
        DataFrame with columns: [group_cols], metric, point_estimate, ci_lower, ci_upper, n
    """
    results = []

    grouped = df.groupby(group_cols)

    for group_keys, group_df in tqdm(
        grouped,
        desc=f"Computing bootstrap CIs for {metric_col}",
        disable=len(grouped) < 10
    ):
        values = group_df[metric_col].values

        lower, upper, point = compute_bootstrap_ci(
            values, n_bootstrap=n_bootstrap, ci=ci, seed=seed
        )

        result = {
            metric_col: metric_col,
            'point_estimate': point,
            'ci_lower': lower,
            'ci_upper': upper,
            'n': len(values[~np.isnan(values)]),
            'n_bootstrap': n_bootstrap
        }

        # Handle both single and multiple group columns
        if len(group_cols) == 1:
            result[group_cols[0]] = group_keys
        else:
            for col, key in zip(group_cols, group_keys):
                result[col] = key

        results.append(result)

    return pd.DataFrame(results)


# =============================================================================
# Significance Tests
# =============================================================================

def paired_wilcoxon_test(
    values1: np.ndarray,
    values2: np.ndarray,
    alternative: str = 'two-sided'
) -> Tuple[float, float, float]:
    """
    Paired Wilcoxon signed-rank test between two sets of values.

    Args:
        values1: First set of values (e.g., tool1 F1 scores)
        values2: Second set of values (e.g., tool2 F1 scores)
        alternative: 'two-sided', 'less', or 'greater'

    Returns:
        Tuple of (statistic, p_value, effect_size_r)
    """
    # Remove NaN pairs
    valid_mask = ~(np.isnan(values1) | np.isnan(values2))
    v1 = values1[valid_mask]
    v2 = values2[valid_mask]

    if len(v1) < 3:
        return np.nan, np.nan, np.nan

    try:
        statistic, p_value = wilcoxon(v1, v2, alternative=alternative)

        # Effect size r = Z / sqrt(N)
        # Approximate Z from statistic for large samples
        n = len(v1)
        z_approx = (statistic - n * (n + 1) / 4) / np.sqrt(n * (n + 1) * (2 * n + 1) / 24)
        effect_size_r = abs(z_approx) / np.sqrt(n)

        return statistic, p_value, effect_size_r

    except ValueError:
        # All differences are zero
        return np.nan, 1.0, 0.0


def permutation_test(
    values1: np.ndarray,
    values2: np.ndarray,
    n_perm: int = 10000,
    seed: Optional[int] = None
) -> Tuple[float, np.ndarray]:
    """
    Permutation test for robust significance testing.

    Args:
        values1: First set of values
        values2: Second set of values
        n_perm: Number of permutations
        seed: Random seed

    Returns:
        Tuple of (p_value, null_distribution)
    """
    if seed is not None:
        np.random.seed(seed)

    # Remove NaN values
    v1 = values1[~np.isnan(values1)]
    v2 = values2[~np.isnan(values2)]

    if len(v1) == 0 or len(v2) == 0:
        return np.nan, np.array([])

    observed_diff = np.mean(v1) - np.mean(v2)

    # Pool and permute
    pooled = np.concatenate([v1, v2])
    n1 = len(v1)

    null_dist = []
    for _ in range(n_perm):
        np.random.shuffle(pooled)
        perm_diff = np.mean(pooled[:n1]) - np.mean(pooled[n1:])
        null_dist.append(perm_diff)

    null_dist = np.array(null_dist)

    # Two-sided p-value
    p_value = np.mean(np.abs(null_dist) >= np.abs(observed_diff))

    return p_value, null_dist


def compute_effect_sizes(
    values1: np.ndarray,
    values2: np.ndarray
) -> Dict[str, float]:
    """
    Compute effect sizes between two groups.

    Args:
        values1: First set of values
        values2: Second set of values

    Returns:
        Dictionary with cohens_d, cohens_d_ci_lower, cohens_d_ci_upper, cliff_delta
    """
    v1 = values1[~np.isnan(values1)]
    v2 = values2[~np.isnan(values2)]

    if len(v1) < 2 or len(v2) < 2:
        return {
            'cohens_d': np.nan,
            'cohens_d_ci_lower': np.nan,
            'cohens_d_ci_upper': np.nan,
            'cliff_delta': np.nan
        }

    # Cohen's d
    mean1, mean2 = np.mean(v1), np.mean(v2)
    var1, var2 = np.var(v1, ddof=1), np.var(v2, ddof=1)
    n1, n2 = len(v1), len(v2)

    # Pooled standard deviation
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

    if pooled_std > 0:
        cohens_d = (mean1 - mean2) / pooled_std
    else:
        cohens_d = 0.0

    # 95% CI for Cohen's d (approximate)
    # SE of d ≈ sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
    se_d = np.sqrt((n1 + n2) / (n1 * n2) + cohens_d**2 / (2 * (n1 + n2)))
    z_95 = 1.96
    cohens_d_ci_lower = cohens_d - z_95 * se_d
    cohens_d_ci_upper = cohens_d + z_95 * se_d

    # Cliff's delta (non-parametric effect size)
    # delta = #(x1 > x2) - #(x1 < x2) / (n1 * n2)
    greater = 0
    less = 0
    for x1 in v1:
        greater += np.sum(x1 > v2)
        less += np.sum(x1 < v2)

    cliff_delta = (greater - less) / (n1 * n2)

    return {
        'cohens_d': cohens_d,
        'cohens_d_ci_lower': cohens_d_ci_lower,
        'cohens_d_ci_upper': cohens_d_ci_upper,
        'cliff_delta': cliff_delta
    }


def compute_all_pairwise_tests(
    df: pd.DataFrame,
    metric_col: str,
    tool_col: str = 'tool',
    group_col: str = 'modification_type',
    n_perm: int = 10000,
    seed: int = 42
) -> pd.DataFrame:
    """
    Compute pairwise significance tests between all tools.

    Args:
        df: DataFrame with metrics
        metric_col: Column name for the metric
        tool_col: Column name for tool identifier
        group_col: Column to stratify by (e.g., modification type)
        n_perm: Number of permutations
        seed: Random seed

    Returns:
        DataFrame with pairwise test results
    """
    results = []
    tools = sorted(df[tool_col].unique())
    groups = df[group_col].unique()

    for group in groups:
        group_df = df[df[group_col] == group]

        for i, tool1 in enumerate(tools):
            for tool2 in tools[i+1:]:
                v1 = group_df[group_df[tool_col] == tool1][metric_col].values
                v2 = group_df[group_df[tool_col] == tool2][metric_col].values

                # Wilcoxon test (paired)
                stat_w, p_w, eff_w = paired_wilcoxon_test(v1, v2)

                # Permutation test
                p_perm, _ = permutation_test(v1, v2, n_perm=n_perm, seed=seed)

                # Effect sizes
                effects = compute_effect_sizes(v1, v2)

                results.append({
                    'tool1': tool1,
                    'tool2': tool2,
                    group_col: group,
                    'metric': metric_col,
                    'test_type': 'wilcoxon',
                    'statistic': stat_w,
                    'p_value': p_w,
                    'effect_size_r': eff_w,
                    'cohens_d': effects['cohens_d'],
                    'cohens_d_ci_lower': effects['cohens_d_ci_lower'],
                    'cohens_d_ci_upper': effects['cohens_d_ci_upper'],
                    'cliff_delta': effects['cliff_delta']
                })

                results.append({
                    'tool1': tool1,
                    'tool2': tool2,
                    group_col: group,
                    'metric': metric_col,
                    'test_type': 'permutation',
                    'statistic': np.nan,
                    'p_value': p_perm,
                    'effect_size_r': effects['cliff_delta'],
                    'cohens_d': effects['cohens_d'],
                    'cohens_d_ci_lower': effects['cohens_d_ci_lower'],
                    'cohens_d_ci_upper': effects['cohens_d_ci_upper'],
                    'cliff_delta': effects['cliff_delta']
                })

    return pd.DataFrame(results)


# =============================================================================
# FDR Correction
# =============================================================================

def apply_fdr_correction(
    p_values: np.ndarray,
    method: str = 'benjamini-hochberg',
    alpha: float = 0.05
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Apply FDR correction to p-values using native implementation.

    Args:
        p_values: Array of p-values
        method: Correction method ('benjamini-hochberg' supported)
        alpha: Significance level

    Returns:
        Tuple of (adjusted_p_values, rejected, significance_level)
    """
    # Use native BH implementation
    rejected, adj_p = benjamini_hochberg_correction(p_values, alpha)
    return adj_p, rejected, adj_p


# =============================================================================
# Main Analysis Functions
# =============================================================================

def run_bootstrap_analysis(
    summary_df: pd.DataFrame,
    by_comparison_df: pd.DataFrame,
    metrics: List[str],
    n_bootstrap: int = 1000,
    ci: float = 95
) -> pd.DataFrame:
    """
    Run bootstrap CI analysis for all metrics.

    Args:
        summary_df: Accuracy summary DataFrame
        by_comparison_df: Per-comparison breakdown DataFrame
        metrics: List of metric columns to analyze
        n_bootstrap: Number of bootstrap iterations
        ci: Confidence interval width

    Returns:
        DataFrame with bootstrap CIs for all metrics
    """
    all_results = []

    # Analysis by tool and modification_type
    for metric in metrics:
        if metric not in summary_df.columns:
            logger.warning(f"Metric {metric} not found in summary, skipping")
            continue

        # Overall by tool
        logger.info(f"Computing bootstrap CIs for {metric} (by tool)")
        result = compute_bootstrap_ci_grouped(
            summary_df,
            metric,
            group_cols=['tool'],
            n_bootstrap=n_bootstrap,
            ci=ci
        )
        result['modification_type'] = 'overall'
        all_results.append(result)

        # By tool and modification_type
        if 'modification_type' in summary_df.columns:
            logger.info(f"Computing bootstrap CIs for {metric} (by tool and modification)")
            result = compute_bootstrap_ci_grouped(
                summary_df,
                metric,
                group_cols=['tool', 'modification_type'],
                n_bootstrap=n_bootstrap,
                ci=ci
            )
            all_results.append(result)

    return pd.concat(all_results, ignore_index=True)


def run_significance_analysis(
    by_comparison_df: pd.DataFrame,
    metrics: List[str],
    n_perm: int = 10000,
    fdr_method: str = 'benjamini-hochberg',
    alpha: float = 0.05
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run pairwise significance tests with FDR correction.

    Args:
        by_comparison_df: Per-comparison metrics DataFrame
        metrics: List of metric columns to test
        n_perm: Number of permutations
        fdr_method: FDR correction method
        alpha: Significance level

    Returns:
        Tuple of (significance_tests_df, fdr_corrected_df)
    """
    all_tests = []
    group_col = 'modification_type' if 'modification_type' in by_comparison_df.columns else 'comparison'

    for metric in metrics:
        if metric not in by_comparison_df.columns:
            logger.warning(f"Metric {metric} not found in comparison data, skipping")
            continue

        logger.info(f"Computing pairwise tests for {metric}")
        tests_df = compute_all_pairwise_tests(
            by_comparison_df,
            metric,
            tool_col='tool',
            group_col=group_col,
            n_perm=n_perm
        )
        all_tests.append(tests_df)

    if not all_tests:
        return pd.DataFrame(), pd.DataFrame()

    all_tests_df = pd.concat(all_tests, ignore_index=True)

    # Apply FDR correction
    logger.info(f"Applying {fdr_method} FDR correction")
    adj_p, rejected, _ = apply_fdr_correction(
        all_tests_df['p_value'].values,
        method=fdr_method,
        alpha=alpha
    )

    all_tests_df['adj_p_value'] = adj_p
    all_tests_df['significant'] = rejected

    # Create FDR summary (use the actual group column name present in the data)
    fdr_summary = all_tests_df[['tool1', 'tool2', group_col, 'metric',
                                 'test_type', 'p_value', 'adj_p_value', 'significant']].copy()

    return all_tests_df, fdr_summary


def run_effect_size_analysis(
    by_comparison_df: pd.DataFrame,
    metrics: List[str]
) -> pd.DataFrame:
    """
    Compute effect sizes between all tool pairs.

    Args:
        by_comparison_df: Per-comparison metrics DataFrame
        metrics: List of metric columns

    Returns:
        DataFrame with effect sizes
    """
    all_effects = []

    group_col = 'modification_type' if 'modification_type' in by_comparison_df.columns else 'comparison'
    tools = sorted(by_comparison_df['tool'].unique())
    groups = by_comparison_df[group_col].unique()

    for metric in metrics:
        if metric not in by_comparison_df.columns:
            continue

        for group in groups:
            group_df = by_comparison_df[by_comparison_df[group_col] == group]

            for i, tool1 in enumerate(tools):
                for tool2 in tools[i+1:]:
                    v1 = group_df[group_df['tool'] == tool1][metric].values
                    v2 = group_df[group_df['tool'] == tool2][metric].values

                    effects = compute_effect_sizes(v1, v2)

                    all_effects.append({
                        'tool1': tool1,
                        'tool2': tool2,
                        group_col: group,
                        'metric': metric,
                        **effects
                    })

    return pd.DataFrame(all_effects)


# =============================================================================
# Main Entry Point
# =============================================================================

def main():
    """Main entry point for benchmark statistics analysis."""

    # Parse Snakemake parameters if available
    try:
        # Snakemake input/output
        summary_input = snakemake.input.summary
        by_comparison_input = snakemake.input.by_comparison
        output_ci = snakemake.output.ci
        output_significance = snakemake.output.significance
        output_fdr = snakemake.output.fdr
        output_effects = snakemake.output.effects

        n_bootstrap = snakemake.params.get('n_bootstrap', 1000)
        alpha = snakemake.params.get('alpha', 0.05)
        fdr_method = snakemake.params.get('fdr_method', 'benjamini-hochberg')

        logger.info("Running in Snakemake mode")

    except NameError:
        # Command line mode
        import argparse

        parser = argparse.ArgumentParser(description='Benchmark Statistical Analysis')
        parser.add_argument('--summary', required=True, help='Accuracy summary TSV')
        parser.add_argument('--by-comparison', required=True, help='Per-comparison TSV')
        parser.add_argument('--output-ci', required=True, help='Output bootstrap CI TSV')
        parser.add_argument('--output-significance', required=True, help='Output significance TSV')
        parser.add_argument('--output-fdr', required=True, help='Output FDR TSV')
        parser.add_argument('--output-effects', required=True, help='Output effect sizes TSV')
        parser.add_argument('--n-bootstrap', type=int, default=1000, help='Bootstrap iterations')
        parser.add_argument('--alpha', type=float, default=0.05, help='Significance level')
        parser.add_argument('--fdr-method', default='benjamini-hochberg', help='FDR method')

        args = parser.parse_args()

        summary_input = args.summary
        by_comparison_input = args.by_comparison
        output_ci = args.output_ci
        output_significance = args.output_significance
        output_fdr = args.output_fdr
        output_effects = args.output_effects
        n_bootstrap = args.n_bootstrap
        alpha = args.alpha
        fdr_method = args.fdr_method

        logger.info("Running in command line mode")

    # Create output directory
    for output_path in [output_ci, output_significance, output_fdr, output_effects]:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Load data
    logger.info(f"Loading summary from: {summary_input}")
    summary_df = pd.read_csv(summary_input, sep='\t')

    logger.info(f"Loading comparison data from: {by_comparison_input}")
    by_comparison_df = pd.read_csv(by_comparison_input, sep='\t')

    # Define metrics to analyze
    metrics = ['precision', 'recall', 'f1', 'auprc', 'auroc', 'mcc', 'specificity']
    metrics = [m for m in metrics if m in summary_df.columns]

    logger.info(f"Analyzing metrics: {metrics}")

    # Run bootstrap analysis
    logger.info("=" * 60)
    logger.info("PHASE 1: Bootstrap Confidence Intervals")
    logger.info("=" * 60)
    bootstrap_df = run_bootstrap_analysis(
        summary_df, by_comparison_df, metrics,
        n_bootstrap=n_bootstrap, ci=95
    )

    # Reorder columns for output
    ci_cols = ['tool', 'modification_type', 'metric', 'point_estimate',
               'ci_lower', 'ci_upper', 'n', 'n_bootstrap']
    ci_cols = [c for c in ci_cols if c in bootstrap_df.columns]
    bootstrap_df = bootstrap_df[ci_cols]
    bootstrap_df.to_csv(output_ci, sep='\t', index=False)
    logger.info(f"Saved bootstrap CIs to: {output_ci}")

    # Run significance analysis
    logger.info("=" * 60)
    logger.info("PHASE 2: Significance Tests")
    logger.info("=" * 60)
    significance_df, fdr_df = run_significance_analysis(
        by_comparison_df, metrics,
        n_perm=10000, fdr_method=fdr_method, alpha=alpha
    )

    if not significance_df.empty:
        significance_df.to_csv(output_significance, sep='\t', index=False)
        logger.info(f"Saved significance tests to: {output_significance}")

        fdr_df.to_csv(output_fdr, sep='\t', index=False)
        logger.info(f"Saved FDR correction to: {output_fdr}")

    # Run effect size analysis
    logger.info("=" * 60)
    logger.info("PHASE 3: Effect Sizes")
    logger.info("=" * 60)
    effects_df = run_effect_size_analysis(by_comparison_df, metrics)

    if not effects_df.empty:
        effects_df.to_csv(output_effects, sep='\t', index=False)
        logger.info(f"Saved effect sizes to: {output_effects}")

    # Summary statistics
    logger.info("=" * 60)
    logger.info("SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Bootstrap CI records: {len(bootstrap_df)}")
    logger.info(f"Significance test records: {len(significance_df)}")
    logger.info(f"Effect size records: {len(effects_df)}")

    if not significance_df.empty:
        n_significant = significance_df['significant'].sum()
        n_total = len(significance_df)
        logger.info(f"Significant comparisons: {n_significant}/{n_total} ({100*n_significant/n_total:.1f}%)")

    logger.info("Statistical analysis complete!")


if __name__ == '__main__':
    main()
