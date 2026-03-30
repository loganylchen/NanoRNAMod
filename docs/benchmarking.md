# Benchmarking System

NanoRNAMod includes a benchmarking system that evaluates RNA modification detection tools against a known truth set. It is activated by setting `benchmark.truth_set` in `config/config.yaml`.

## Output Directory Structure

All outputs are under `{project}/results/benchmarks/`.

```
benchmarks/
├── cross_tool/                          # Aggregated cross-tool comparison
│   ├── accuracy_summary.tsv             # All tool × comparison, best score column
│   ├── accuracy_summary_overall.tsv     # Per-tool aggregated metrics (mean ± std)
│   ├── accuracy_summary_by_comparison.tsv  # Per-tool × comparison × mode (native + fair)
│   ├── by_tool.tsv                      # Final per-tool summary with coverage data
│   └── best_scores.tsv                  # Score column selection rationale per tool
│
├── per_tool/{tool}/                     # Per-tool results
│   ├── {comparison}/
│   │   ├── native/                      # Evaluation on tool's own reported sites
│   │   │   ├── score_comparison.tsv     # All score columns evaluated
│   │   │   ├── best_metrics.tsv         # Best score column metrics
│   │   │   └── best_score.txt           # Name of best score column
│   │   ├── fair/                        # Evaluation on union of ALL tools' sites
│   │   │   ├── score_comparison.tsv
│   │   │   ├── best_metrics.tsv
│   │   │   └── best_score.txt
│   │   └── figures/                     # Per-tool × comparison figures
│   │       ├── roc_curve.pdf            # ROC curve (native mode)
│   │       ├── pr_curve.pdf             # Precision-recall curve
│   │       ├── score_distribution.pdf   # Score histogram by label (positive/negative)
│   │       ├── native_vs_fair.pdf       # Bar chart comparing native vs fair metrics
│   │       └── *_data.tsv              # Source data for each figure
│   ├── {sample}/                        # Per-sample tools (tandemmod, etc.)
│   │   ├── native/                      # Native evaluation only (no fair mode)
│   │   └── figures/
│   └── figures/                         # Per-tool summary (across all comparisons)
│       ├── accuracy_by_comparison.pdf   # AUROC/F1 grouped by comparison
│       ├── native_vs_fair_summary.pdf   # Paired native vs fair per comparison
│       ├── score_columns_comparison.pdf # All score columns ranked by AUROC
│       ├── threshold_curve.pdf          # F1 vs threshold for best score column
│       └── *_data.tsv
│
├── figures/                             # Publication figures (flat layout)
│   ├── fig1_overall_accuracy.pdf        # F1 score comparison with bootstrap CI
│   ├── fig2_precision_recall.pdf        # Precision vs recall trade-off
│   ├── fig3_auroc.pdf                   # AUROC comparison with CI
│   ├── fig4_native_vs_fair.pdf          # Native vs fair AUROC per tool
│   ├── fig5_best_score.pdf              # Best score column selection dot plot
│   ├── fig6_coverage_sensitivity.pdf    # Performance vs coverage depth
│   ├── fig7_resource_usage.pdf          # Runtime vs memory (bubble = F1)
│   ├── fig8_tool_ranking.pdf            # Lollipop: tools ranked by AUROC
│   ├── sfig1_per_comparison.pdf         # F1 by comparison (supplementary)
│   ├── sfig2_native_vs_fair_detail.pdf  # Per-tool native/fair breakdown
│   ├── sfig3_threshold_robustness.pdf   # Threshold stability (CV)
│   ├── sfig4_effect_sizes.pdf           # Pairwise Cohen's d
│   ├── sfig5_score_heatmap.pdf          # Score column AUROC heatmap
│   └── fig*_data.tsv                    # Co-located source data for each figure
│
├── coverage/                            # Coverage analysis (intermediate)
│   └── {comparison}/
│       ├── {tool}_covered.tsv           # Truth sites covered by each tool
│       ├── union.tsv                    # Union of truth-covered sites
│       ├── union_predictions.tsv        # Union of ALL tool-reported sites
│       └── called_sites.tsv             # Per-tool site count
│
├── statistics/                          # Statistical analysis (intermediate)
│   ├── bootstrap_ci.tsv                 # 95% bootstrap confidence intervals
│   ├── significance_tests.tsv           # Paired Wilcoxon tests
│   ├── fdr_corrected.tsv                # FDR-corrected p-values
│   ├── effect_sizes.tsv                 # Cohen's d, Cliff's delta
│   └── auc_comparison.tsv              # Pairwise AUC comparison (bootstrap DeLong-like test)
│
├── sensitivity/                         # Sensitivity analysis (intermediate)
│   ├── coverage_analysis.tsv            # Performance stratified by read depth
│   ├── score_distribution.tsv           # Score separation metrics
│   └── threshold_robustness.tsv         # Threshold stability across CV splits
│
├── benchmark_report.pdf                 # Comprehensive PDF report
├── resource_summary.tsv                 # Per-rule resource usage
├── resource_by_tool.tsv                 # Resource usage grouped by tool
└── resource_by_tool.pdf                 # Resource usage figure
```

## Key Concepts

### Native vs Fair Evaluation

Each comparison tool is evaluated in two modes:

- **Native mode**: Evaluates the tool on its own reported sites. Only positions the tool actually reports are scored. This shows the tool's accuracy on sites it is confident enough to call.

- **Fair mode**: Evaluates the tool on the union of ALL positions reported by ANY tool (`union_predictions.tsv`). Sites a tool does not report receive the worst possible score (p-values → 1.0, regular scores → minimum - 1). This penalizes tools that miss sites and enables apples-to-apples comparison across tools.

### Score Column Selection

Each tool may output multiple score columns (e.g., p-value, effect size, z-score). The benchmarking system evaluates all candidate columns and selects the best one per tool using:

```
selection_score = mean_AUROC - std_AUROC
```

This penalizes inconsistent columns — a column with AUROC 0.75 ± 0.01 beats 0.76 ± 0.15.

### P-value Handling

- Raw p-values (lower = more significant) are transformed to `-log10(p)` so that higher = better, matching the convention for other score types.
- Some tools (e.g., differr) output pre-transformed `-log10 P value` columns — these are used directly without additional transformation.
- `inf` values (from `-log10(0)`) are capped at 1000.0.
- In fair mode, uncalled p-value sites are filled with 1.0 (not significant), so `-log10(1) = 0`.

### Metric Definitions

| Metric | Description |
|--------|-------------|
| **AUROC** | Area Under ROC Curve — threshold-independent discrimination ability (0.5 = random, 1.0 = perfect) |
| **PRAUC** | Area Under Precision-Recall Curve — useful when classes are imbalanced |
| **F1** | Harmonic mean of precision and recall at the optimal threshold |
| **Precision** | Fraction of positive predictions that are correct (TP / (TP + FP)) |
| **Recall** | Fraction of true positives that are detected (TP / (TP + FN)) |
| **Youden's J** | Sensitivity + Specificity - 1 (optimal point on ROC curve, range -1 to 1) |

### Threshold Selection

Each score column is evaluated with two threshold selection methods, computed independently for native and fair modes:

| Method | Criterion | Use case |
|--------|-----------|----------|
| **F1-optimal** (`best_threshold`) | Maximizes F1 = 2PR/(P+R) | Balanced precision-recall trade-off |
| **Youden-optimal** (`youden_threshold`) | Maximizes TPR - FPR (Youden's J) | Optimal point on the ROC curve |

The `_original` variants give the threshold in the original scale (before `-log10` transform for p-values).

**Native vs Fair thresholds**: These are computed separately because the evaluation sites differ:
- **Native**: Only sites the tool reports (smaller n, high-confidence sites)
- **Fair**: Union of all tools' reported sites (larger n, includes uncalled sites filled with worst-case scores)

Comparing thresholds across modes helps decide whether to use a tool's native output or to apply it within a multi-tool ensemble.

### AUC Comparison (Bootstrap DeLong-like Test)

The `statistics/auc_comparison.tsv` file contains pairwise comparisons of AUROC and PRAUC between all tool pairs. For each pair:

- **observed_diff**: Mean AUC difference (tool1 - tool2) across comparisons
- **ci_lower / ci_upper**: 95% bootstrap confidence interval on the difference
- **p_value**: Two-sided bootstrap p-value testing H0: no difference
- **adj_p_value**: FDR-corrected p-value (Benjamini-Hochberg)
- **significant**: Whether the difference is significant after FDR correction

This is analogous to the DeLong test but operates on aggregated AUC values across comparisons using paired bootstrap resampling, rather than on raw prediction scores.

## Primary Output Files

These are the files most useful for analysis:

| File | Purpose |
|------|---------|
| `cross_tool/by_tool.tsv` | **Start here.** Per-tool summary: best score column, AUROC, F1, precision, recall, coverage |
| `cross_tool/accuracy_summary_by_comparison.tsv` | Per-tool × comparison metrics with native/fair mode breakdown |
| `cross_tool/best_scores.tsv` | Why each score column was selected (all candidates with mean AUROC) |
| `figures/fig2_precision_recall.pdf` | Precision, recall, and PRAUC comparison |
| `figures/fig4_native_vs_fair.pdf` | Visual comparison of native vs fair AUROC |
| `figures/fig8_tool_ranking.pdf` | Tools ranked by overall performance |
| `per_tool/{tool}/figures/` | Detailed per-tool analysis across comparisons |
| `statistics/auc_comparison.tsv` | Pairwise AUC significance tests between tools |
| `benchmark_report.pdf` | Comprehensive multi-page report |
