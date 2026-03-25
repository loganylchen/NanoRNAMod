# NanoRNAMod Benchmarking Module Design Document

## 1. Overview

### 1.1 Purpose
The benchmarking module provides standardized accuracy and performance evaluation for RNA modification detection tools integrated into NanoRNAMod. It enables:
- **Accuracy benchmarking**: Compare tool predictions against known ground truth modification sites
- **Resource benchmarking**: Track computational resources (CPU, memory, I/O) consumed by each tool
- **Cross-tool comparison**: Generate unified reports comparing all activated tools

### 1.2 Current Implementation Status

The following components are implemented:

**Core Accuracy Rules & Scripts:**
- `workflow/rules/benchmark_accuracy.smk` - Accuracy evaluation rules
- `workflow/rules/benchmark_report.smk` - Resource aggregation rule
- `workflow/rules/benchmark_viz.smk` - Publication-quality visualization
- `workflow/scripts/accuracy_benchmark.py` - Core accuracy metrics computation
- `workflow/scripts/aggregate_benchmarks.py` - Resource benchmark aggregation
- `workflow/scripts/benchmark_utils.py` - Shared utility functions

**Statistical Analysis:**
- `workflow/scripts/benchmark_statistics.py` - Bootstrap CIs, Wilcoxon tests, FDR correction

**Sensitivity Analysis:**
- `workflow/scripts/benchmark_sensitivity.py` - Coverage-stratified analysis, cross-validation

**Threshold Optimization:**
- `workflow/scripts/benchmark_threshold.py` - Optimal threshold finding per tool
- `workflow/scripts/benchmark_score_optimization.py` - Multi-threshold evaluation
- `workflow/scripts/benchmark_multithreshold.py` - Threshold sensitivity analysis

**Detailed Reports:**
- `workflow/scripts/benchmark_detailed.py` - Site-by-site prediction analysis

**Visualization:**
- `workflow/scripts/benchmark_plots.py` - Python visualization (deprecated, use R)
- `workflow/scripts/benchmark_pdf_report.py` - PDF report generation
- `workflow/scripts/benchmark_resource_by_tool.py` - Resource usage plots
- R scripts in `workflow/scripts/` for ggplot2-based publication figures

---

## 2. Input Data Formats

### 2.1 Truth Set Format (TSV)

**File path**: Configured via `benchmark.truth_set` in `config/config.yaml`

**Required columns**:
| Column | Type | Description |
|--------|------|-------------|
| `transcript` | string | Transcript identifier (must match tool output) |
| `position` | integer | 0-based position on transcript |
| `modification_type` | string | Short code (e.g., `m6A`, `Psi`, `m5C`) |
| `modification_name` | string | Full name (e.g., `N6-methyladenosine`) - optional |
| `label` | string | `positive` for true sites, `-` for negative sites (optional) |

**Example**:
```tsv
transcript	position	modification_type	modification_name	label
ecoli16S	516	Psi	pseudouridine	positive
ecoli16S	527	m7G	N7-methylguanosine	positive
ecoli16S	1000	m6A	N6-methyladenosine	-
```

**Label semantics**:
- `positive` (or any non-`-` value): Confirmed modification site
- `-`: Confirmed negative site (no modification at this position)
- Missing label column: All sites treated as positive; negatives auto-inferred from predictions

### 2.2 Tool Output Format (Standardized)

All tools must output a TSV file at:
```
{project}/results/modifications/{tool}/{sample_or_comparison}/{tool}_results.tsv
```

**Required columns** (normalized by post-processing scripts):
| Column | Aliases | Description |
|--------|---------|-------------|
| `transcript` | `id`, `ref_id`, `chrom` | Transcript identifier |
| `position` | `pos`, `start`, `start_loc`, `transcript_pos` | 0-based position |

**Optional score columns** (used for ranking metrics):
| Tool | Primary Score Column | Alternative Columns |
|------|---------------------|---------------------|
| xpore | `p_value` | `diff_mod`, `diff_mod_frac` |
| nanocompore | `pvalue` | `logit_pvalue`, `logit` |
| baleen | `mod_score` | `score`, `kmer_score` |
| differr | `-log10 P value` | `score`, `pvalue` |
| drummer | `p_value` | `z_score` |
| eligos2 | `pvalue` | `esb`, `oddsR` |
| epinano | `z_score_prediction` | `delta_sum_err` |
| psipore | `pvalue` | `p_value`, `score` |
| tandemmod | `probability` | `prob`, `score` |
| directrm | `probability` | `prob`, `mod_prob` |
| m6atm | `stoichiometry` | `probability`, `prob` |
| rnano | `probability` | `score`, `prob` |
| nanopsu | `pvalue` | `p_value`, `score` |
| nanomud | `probability` | `pvalue`, `score` |
| penguin | `probability` | `score`, `pvalue` |
| pybaleen | `mod_score` | `score` |

---

## 3. Evaluation Metrics

### 3.1 Binary Classification Metrics

Metrics are computed at multiple positional tolerance windows (configured via `benchmark.window`).

| Metric | Formula | Description |
|--------|---------|-------------|
| **Precision** | TP / (TP + FP) | Fraction of predictions that are correct |
| **Recall (Sensitivity)** | TP / (TP + FN) | Fraction of truth sites detected |
| **F1 Score** | 2 × (P × R) / (P + R) | Harmonic mean of precision and recall |
| **Specificity** | TN / (TN + FP) | Fraction of negative sites correctly rejected |
| **MCC** | See formula below | Matthews Correlation Coefficient |

**MCC Formula**:
```
MCC = (TP × TN - FP × FN) / √((TP+FP)(TP+FN)(TN+FP)(TN+FN))
```

**Confusion Matrix Definitions**:
- **TP**: Truth positive site with prediction within ±window nucleotides
- **FN**: Truth positive site with no prediction within ±window
- **FP**: Prediction with no truth site within ±window
- **TN**: Truth negative site with no prediction within ±window

### 3.2 Ranking Metrics (Score-based)

| Metric | Description | Interpretation |
|--------|-------------|----------------|
| **AUROC** | Area Under ROC Curve | Probability that a random positive ranks higher than random negative |
| **AUPRC** | Area Under Precision-Recall Curve | Better for imbalanced datasets |

**Score Transformation**:
- P-value columns: Transformed to `-log10(p_value)` (higher = better)
- Probability/score columns: Used directly

### 3.3 Summary Statistics

| Metric | Description |
|--------|-------------|
| `called_sites` | Total predictions made by tool |
| `total_truth` | Total positive truth sites |
| `total_predicted` | Same as `called_sites` |
| `total_negative` | Total negative truth sites (if provided) |

### 3.4 Statistical Analysis (Publication Rigor)

#### 3.4.1 Bootstrap Confidence Intervals

Non-parametric bootstrap resampling provides confidence intervals for all metrics:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_bootstrap` | 1000 | Number of bootstrap iterations |
| `alpha` | 0.05 | Significance level (95% CI) |

**Method**: Percentile-based CIs from bootstrap distribution of each metric.

#### 3.4.2 Significance Testing

**Wilcoxon Rank-Sum Tests**: Compare each tool against a reference baseline
- Tests for significant differences in precision, recall, F1
- Multiple testing correction via FDR

| FDR Method | Description |
|------------|-------------|
| `benjamini-hochberg` | Standard FDR control (default) |
| `bonferroni` | Conservative family-wise error rate |
| `holm` | Step-down Bonferroni |

#### 3.4.3 Output Files

| File | Content |
|------|---------|
| `statistics/summary_statistics.tsv` | Bootstrap CIs for all metrics |
| `statistics/pairwise_tests.tsv` | Significance tests with FDR correction |
| `statistics/statistics_report.txt` | Human-readable summary |

### 3.5 Sensitivity Analysis

#### 3.5.1 Coverage-Stratified Analysis

Performance evaluation stratified by read coverage depth:

```yaml
coverage_bins: [0, 10, 20, 50, 100, 200, 500]
```

For each bin, computes precision, recall, F1 to identify:
- Coverage thresholds for reliable detection
- Performance degradation at low coverage
- Tool-specific coverage dependencies

#### 3.5.2 Threshold Robustness (Cross-Validation)

K-fold cross-validation for threshold stability:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_splits` | 5 | Number of CV folds |

Evaluates how consistently each tool's optimal threshold generalizes across data splits.

#### 3.5.3 Output Files

| File | Content |
|------|---------|
| `sensitivity/coverage_stratified.tsv` | Metrics by coverage bin |
| `sensitivity/cross_validation.tsv` | Threshold robustness across folds |
| `sensitivity/sensitivity_report.txt` | Analysis summary |

### 3.6 Publication Figures

#### 3.6.1 Main Figures

| Figure | Description | Format |
|--------|-------------|--------|
| `fig1_overall_accuracy.pdf` | Bar plot of F1 scores with CIs | PDF, Nature theme |
| `fig2_pr_curves.pdf` | Precision-recall curves per tool | PDF, Nature theme |
| `fig3_roc_curves.pdf` | ROC curves with AUROC | PDF, Nature theme |
| `fig4_window_sensitivity.pdf` | F1 vs positional tolerance | PDF, Nature theme |
| `fig5_tool_comparison.pdf` | Pairwise comparison heatmap | PDF, Nature theme |
| `fig6_resource_usage.pdf` | CPU/memory/runtime comparison | PDF, Nature theme |

#### 3.6.2 Supplementary Figures

| Figure | Description | Format |
|--------|-------------|--------|
| `fig_s1_score_distributions.pdf` | Score distributions per tool | PDF |
| `fig_s2_mod_type_breakdown.pdf` | Metrics by modification type | PDF |
| `fig_s3_coverage_sensitivity.pdf` | Coverage-stratified analysis | PDF |
| `fig_s4_threshold_robustness.pdf` | Cross-validation stability | PDF |

#### 3.6.3 Source Data Files

All plotting data preserved for reproducibility:

| File | Content |
|------|---------|
| `data/fig1_source_data.tsv` | Data for Figure 1 |
| `data/fig2_source_data.tsv` | Data for Figure 2 |
| ... | (one file per figure) |

#### 3.6.4 Visualization Standards

- **Theme**: ggplot2 with Nature/Cell/Science journal presets
- **Colors**: Okabe-Ito colorblind-friendly palette
- **Fonts**: Helvetica/Arial family, 7-10pt base size
- **Resolution**: 300 DPI minimum for all figures
- **Dimensions**: Single column (3.5"), double column (7.0")

---

## 4. Output File Structure

### 4.1 Directory Layout

```
{project}/results/benchmarks/
├── accuracy_summary.tsv           # Per-modification_type metrics
├── accuracy_summary_overall.tsv   # Aggregated across all mod types
├── resource_summary.tsv           # CPU/memory/IO benchmarks
│
├── statistics/                    # Statistical analysis outputs
│   ├── bootstrap_ci.tsv           # Bootstrap confidence intervals
│   ├── significance_tests.tsv     # Pairwise Wilcoxon tests
│   └── ranking_summary.tsv        # Tool rankings with CIs
│
├── sensitivity/                   # Sensitivity analysis outputs
│   ├── coverage_stratified.tsv    # Metrics by coverage bin
│   ├── cv_stability.tsv           # Cross-validation results
│   └── robustness_summary.tsv     # Threshold robustness analysis
│
├── figures/                       # Publication-quality figures
│   ├── main/                      # Main text figures
│   │   ├── fig1_overall_accuracy.pdf
│   │   ├── fig2_pr_curves.pdf
│   │   ├── fig3_metric_comparison.pdf
│   │   ├── fig4_resource_usage.pdf
│   │   ├── fig5_statistical_comparison.pdf
│   │   └── fig6_sensitivity_analysis.pdf
│   └── supplementary/             # Supplementary figures
│       ├── figS1_*.pdf
│       └── ...
│
├── data/                          # Source data for all figures
│   ├── fig1_source_data.tsv
│   ├── fig2_source_data.tsv
│   └── ...
│
├── optimal_thresholds.tsv         # Best score thresholds per tool
├── threshold_evaluation.tsv       # Metrics at each evaluated threshold
├── score_distributions.tsv        # Score statistics per tool
├── detailed_predictions.tsv       # Site-by-site prediction analysis
└── detailed_truth.tsv             # Site-by-site truth coverage

logs/{project}/
├── accuracy_benchmark/
│   └── accuracy.log
├── aggregate_benchmarks/
│   └── aggregate.log
└── benchmark_viz/
    └── viz.log

benchmarks/{project}/
├── {sample}.{tool}.benchmark.txt  # Per-rule resource usage
└── ...
```

### 4.2 accuracy_summary.tsv Columns

```tsv
tool	modification_type	window	precision	recall	f1	tp	fp	fn	tn	specificity	mcc	auprc	auroc	called_sites	total_truth	total_predicted	total_negative
```

### 4.3 accuracy_summary_overall.tsv Columns

```tsv
tool	window	precision	recall	f1	tp	fp	fn	tn	specificity	mcc	auprc	auroc	called_sites	total_truth	total_predicted	total_negative
```

### 4.4 statistics_summary.tsv Columns (NEW)

Statistical analysis with confidence intervals and significance tests:

```tsv
tool	modification_type	window	metric	value	ci_lower	ci_upper	bootstrap_n
```

### 4.5 pairwise_tests.tsv Columns (NEW)

Pairwise significance tests between tools:

```tsv
tool_a	tool_b	modification_type	window	metric	p_value	p_value_corrected	significant
```

### 4.6 sensitivity_summary.tsv Columns (NEW)

Coverage-stratified performance analysis:

```tsv
tool	modification_type	window	coverage_bin	n_sites	precision	recall	f1
```

### 4.7 threshold_robustness.tsv Columns (NEW)

Cross-validation threshold stability:

```tsv
tool	modification_type	split	optimal_threshold	f1_at_threshold
```

### 4.8 optimal_thresholds.tsv Columns

```tsv
tool	modification_type	window	criterion	optimal_threshold	optimal_f1	precision	recall	tp	fp	fn
```

### 4.9 detailed_predictions.tsv Columns (NEW)

Site-by-site prediction analysis:

```tsv
tool	transcript	position	score	truth_match	distance_to_truth	matched_modification
```

### 4.10 detailed_truth.tsv Columns (NEW)

Truth site detection coverage:

```tsv
tool	transcript	position	modification_type	detected	nearest_prediction_distance	nearest_prediction_score
```

### 4.11 resource_summary.tsv Columns

Based on Snakemake benchmark output:
```tsv
s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load	cpu_time	tool	sample
```

| Column | Description |
|--------|-------------|
| `s` | Wall-clock time in seconds |
| `h:m:s` | Wall-clock time as HH:MM:SS |
| `max_rss` | Maximum resident set size (MB) |
| `max_vms` | Maximum virtual memory (MB) |
| `io_in` | Total I/O read (MB) |
| `io_out` | Total I/O write (MB) |
| `mean_load` | Mean CPU load (%) |
| `cpu_time` | Total CPU time (seconds) |

---

## 5. Integration with Existing Workflow

### 5.1 Configuration Integration

**config/config.yaml**:
```yaml
benchmark:
  truth_set: "config/benchmark_ecoli_rRNA.tsv"  # Path to truth set (empty to skip)
  window: [0, 1, 2, 5]                          # Positional tolerance windows

  # Multi-threshold analysis
  n_thresholds: 50                              # Number of thresholds to evaluate
  custom_thresholds: []                         # Custom thresholds to evaluate (optional)

  # Statistical analysis (publication rigor)
  n_bootstrap: 1000                             # Bootstrap iterations for CIs
  alpha: 0.05                                   # Significance level for tests
  fdr_method: "benjamini-hochberg"              # FDR correction method

  # Sensitivity analysis
  coverage_bins: [0, 10, 20, 50, 100, 200, 500] # Coverage depth bins
  n_splits: 5                                   # Cross-validation splits

  # Publication figure settings
  figure_theme: "nature"                        # Journal theme (nature/cell/science)
  dpi: 300                                      # Figure resolution
  fig_format: "pdf"                             # Output format (pdf/png/svg)
  include_supplementary: true                   # Generate supplementary figures
  data_export_format: "tsv"                     # Format for source data
```

**workflow/schemas/config.schema.yaml** (already defined):
```yaml
benchmark:
  type: object
  properties:
    truth_set:
      type: string
    window:
      oneOf:
        - type: integer
        - type: array
          items:
            type: integer
```

### 5.2 Workflow DAG Integration

```
┌─────────────────────────────────────────────────────────────────┐
│                     Modification Detection                       │
│  ┌─────────┐ ┌─────────────┐ ┌─────────┐ ┌─────────┐            │
│  │ xpore   │ │ nanocompore │ │ baleen  │ │ ...     │            │
│  └────┬────┘ └──────┬──────┘ └────┬────┘ └────┬────┘            │
│       │             │             │           │                  │
│       └─────────────┴──────┬──────┴───────────┘                  │
│                            ▼                                     │
│              ┌─────────────────────────┐                         │
│              │  {tool}_results.tsv     │                         │
│              │  (per tool/comparison)  │                         │
│              └────────────┬────────────┘                         │
└───────────────────────────┼─────────────────────────────────────┘
                            │
                            ▼
┌───────────────────────────────────────────────────────────────────┐
│                      Benchmarking Module                           │
│                                                                    │
│  ┌─────────────────────────────────────────────────────────────┐  │
│  │                   accuracy_benchmark rule                    │  │
│  │  Input:                                                      │  │
│  │    - All activated tool results                              │  │
│  │    - Truth set TSV                                           │  │
│  │  Output:                                                     │  │
│  │    - accuracy_summary.tsv                                    │  │
│  │    - accuracy_summary_overall.tsv                            │  │
│  └─────────────────────────────────────────────────────────────┘  │
│                                                                    │
│  ┌─────────────────────────────────────────────────────────────┐  │
│  │                  aggregate_benchmarks rule                   │  │
│  │  Input:                                                      │  │
│  │    - benchmarks/{project}/*.benchmark.txt                    │  │
│  │  Output:                                                     │  │
│  │    - resource_summary.tsv                                    │  │
│  └─────────────────────────────────────────────────────────────┘  │
│                                                                    │
└───────────────────────────────────────────────────────────────────┘
```

### 5.3 Rule Dependencies

```
accuracy_benchmark:
  depends_on:
    - All modetect_* rules (via get_all_result_tsvs function)
  triggers_when: benchmark.truth_set is non-empty

aggregate_benchmarks:
  depends_on:
    - All rules with benchmark: directive
  always_runs: Yes (aggregates existing benchmark files)
```

### 5.4 Final Output Registration

In `workflow/rules/common.smk` → `get_final_output()`:
```python
# Benchmarking outputs
final_output += [f"{RESULT_ROOT}/benchmarks/resource_summary.tsv"]
if config.get("benchmark", {}).get("truth_set", ""):
    final_output += [f"{RESULT_ROOT}/benchmarks/accuracy_summary.tsv"]
```

---

## 6. Rules and Scripts Inventory

### 6.1 Existing Rules

| Rule File | Rule Name | Purpose |
|-----------|-----------|---------|
| `workflow/rules/benchmark_accuracy.smk` | `accuracy_benchmark` | Compute accuracy metrics against truth set |
| `workflow/rules/benchmark_accuracy.smk` | `benchmark_detailed` | Site-by-site analysis |
| `workflow/rules/benchmark_accuracy.smk` | `benchmark_threshold` | Optimal threshold finding |
| `workflow/rules/benchmark_accuracy.smk` | `benchmark_score_optimization` | Multi-threshold evaluation |
| `workflow/rules/benchmark_report.smk` | `aggregate_benchmarks` | Aggregate resource benchmarks |
| `workflow/rules/benchmark_viz.smk` | `benchmark_r_figures` | Publication-quality R figures |
| `workflow/rules/benchmark_viz.smk` | `benchmark_pdf_report` | Comprehensive PDF report |

### 6.2 Existing Scripts

**Core Scripts:**
| Script File | Purpose |
|-------------|---------|
| `workflow/scripts/accuracy_benchmark.py` | Core accuracy computation logic |
| `workflow/scripts/aggregate_benchmarks.py` | Resource benchmark aggregation |
| `workflow/scripts/benchmark_utils.py` | Shared utility functions |
| `workflow/scripts/benchmark_detailed.py` | Site-by-site comparison reports |

**Statistical Analysis:**
| Script File | Purpose |
|-------------|---------|
| `workflow/scripts/benchmark_statistics.py` | Bootstrap CIs, Wilcoxon tests, FDR correction |

**Sensitivity Analysis:**
| Script File | Purpose |
|-------------|---------|
| `workflow/scripts/benchmark_sensitivity.py` | Coverage stratification, cross-validation |

**Threshold Analysis:**
| Script File | Purpose |
|-------------|---------|
| `workflow/scripts/benchmark_threshold.py` | Optimal threshold finding per tool |
| `workflow/scripts/benchmark_score_optimization.py` | Multi-threshold evaluation |
| `workflow/scripts/benchmark_multithreshold.py` | Threshold sensitivity analysis |

**Visualization:**
| Script File | Purpose |
|-------------|---------|
| `workflow/scripts/benchmark_plots.py` | Python visualization (deprecated) |
| `workflow/scripts/benchmark_pdf_report.py` | PDF report generation |
| `workflow/scripts/benchmark_resource_by_tool.py` | Resource usage plots |

**Negative Control Strategies:**
| Script File | Purpose |
|-------------|---------|
| `workflow/scripts/benchmark_kmer_negatives.py` | K-mer based negative site generation |
| `workflow/scripts/benchmark_same_base_negatives.py` | Same-base negative control strategy |

### 6.3 Conda Environments

| Environment File | Dependencies |
|------------------|--------------|
| `workflow/envs/pandas.yaml` | pandas, numpy, scikit-learn |
| `workflow/envs/r_viz.yaml` | R, ggplot2, dplyr, tidyr, cowplot |

---

## 7. Implemented Enhancements

### 7.1 Per-Tool Detailed Reports ✅ IMPLEMENTED

**Output files**:
- `{project}/results/benchmarks/detailed_predictions.tsv`
- `{project}/results/benchmarks/detailed_truth.tsv`

**Script**: `workflow/scripts/benchmark_detailed.py`

**Sample output**:
```tsv
tool	transcript	position	score	truth_match	distance_to_truth	matched_modification
xpore	ecoli16S	516	0.001	yes	0	Psi
xpore	ecoli16S	527	0.002	yes	0	m7G
nanocompore	ecoli16S	516	0.003	yes	0	Psi
```

### 7.2 Publication-Quality Visualization ✅ IMPLEMENTED

**Rule**: `benchmark_r_figures` in `workflow/rules/benchmark_viz.smk`

**Output directory**: `{project}/results/benchmarks/figures/`

**Main figures**:
- `main/fig1_overall_accuracy.pdf` - Bar plot with CIs
- `main/fig2_pr_curves.pdf` - Precision-recall curves
- `main/fig3_roc_curves.pdf` - ROC curves with AUROC
- `main/fig4_window_sensitivity.pdf` - Positional tolerance analysis
- `main/fig5_tool_comparison.pdf` - Pairwise comparison heatmap
- `main/fig6_resource_usage.pdf` - Resource benchmarks

**Supplementary figures**:
- `supplementary/fig_s1_score_distributions.pdf`
- `supplementary/fig_s2_mod_type_breakdown.pdf`
- `supplementary/fig_s3_coverage_sensitivity.pdf`
- `supplementary/fig_s4_threshold_robustness.pdf`

**Source data**: All plotting data preserved in `data/` subdirectory.

### 7.3 Threshold Optimization ✅ IMPLEMENTED

**Output files**:
- `{project}/results/benchmarks/optimal_thresholds.tsv`
- `{project}/results/benchmarks/threshold_evaluation.tsv`
- `{project}/results/benchmarks/score_distributions.tsv`

**Script**: `workflow/scripts/benchmark_threshold.py`

**Features**:
- Finds threshold maximizing F1, precision, or recall
- Evaluates at multiple windows
- Reports confusion matrix at optimal threshold

### 7.4 Statistical Analysis ✅ IMPLEMENTED

**Output directory**: `{project}/results/benchmarks/statistics/`

**Features**:
- Bootstrap confidence intervals (default: 1000 iterations)
- Pairwise Wilcoxon rank-sum tests
- FDR correction (Benjamini-Hochberg, Bonferroni, Holm)
- Significance stars for publication tables

**Script**: `workflow/scripts/benchmark_statistics.py`

### 7.5 Sensitivity Analysis ✅ IMPLEMENTED

**Output directory**: `{project}/results/benchmarks/sensitivity/`

**Features**:
- Coverage-stratified performance (configurable bins)
- K-fold cross-validation for threshold robustness
- Modification type breakdown

**Script**: `workflow/scripts/benchmark_sensitivity.py`

---

## 8. Future Enhancements

### 8.1 Multi-Truth-Set Support

**Configuration**:
```yaml
benchmark:
  truth_sets:
    - name: "ecoli_rrna"
      path: "config/benchmark_ecoli_rRNA.tsv"
    - name: "hekar_m6a"
      path: "config/benchmark_HEKAR_m6A.tsv"
```

### 8.2 Motif-Context Analysis

**Feature**: Stratify by sequence context:
- DRACH motif for m6A (D=A/G/U, R=A/G, H=A/C/U)
- Hairpin structures for Ψ
- Custom motifs from user

### 8.3 Interactive Dashboard

**Feature**: Web-based exploration of:
- Interactive ROC/PR curves with hover details
- Per-site truth/prediction visualization
- Download filtered subsets

---

## 9. Implementation Checklist

### Phase 1: Core Enhancements ✅ COMPLETE
- [x] Add `modification_name` column support to truth set parsing
- [x] Implement per-tool detailed reports (`benchmark_detailed.py`)
- [x] Add threshold optimization logic (`benchmark_threshold.py`)
- [x] Improve score column detection with tool-specific mappings
- [x] Create shared utility module (`benchmark_utils.py`)

### Phase 2: Statistical Rigor ✅ COMPLETE
- [x] Bootstrap confidence intervals for all metrics
- [x] Wilcoxon rank-sum significance tests
- [x] FDR correction for multiple testing
- [x] Pairwise tool comparison with statistical output

### Phase 3: Sensitivity Analysis ✅ COMPLETE
- [x] Coverage-stratified performance analysis
- [x] Cross-validation for threshold robustness
- [x] Modification type breakdown

### Phase 4: Visualization ✅ COMPLETE
- [x] Create `benchmark_viz.smk` rule file
- [x] Implement R/ggplot2 visualization with Nature theme
- [x] Generate main figures (6 publication figures)
- [x] Generate supplementary figures (4 additional)
- [x] Preserve all source data for reproducibility

### Phase 5: Documentation ⏳ IN PROGRESS
- [x] Update `config/config.yaml` comments
- [x] Update benchmarking design document
- [ ] Add benchmarking section to main README
- [ ] Create example truth sets for common organisms

### Phase 6: Future Enhancements
- [ ] Multi-truth-set support for batch evaluations
- [ ] Motif-context analysis (e.g., DRACH for m6A)
- [ ] Interactive web-based dashboard
- [ ] Integration with GEO/SRA metadata

---

## 10. API Reference

### 10.1 Helper Functions (common.smk)

```python
def get_all_result_tsvs(wildcards):
    """
    Collect all *_results.tsv paths for activated tools.
    Returns list of file paths for both per-comparison and per-sample tools.
    """
```

### 10.2 Shared Utilities (benchmark_utils.py)

```python
def tool_from_path(path: str) -> str:
    """Extract tool name from result file path."""

def comparison_from_path(path: str) -> str:
    """Extract comparison name from result file path."""

def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Map various column names to standard 'transcript' and 'position'."""

def match_column(df_columns: list, pattern: str) -> str:
    """Match a pattern against DataFrame columns using multiple strategies."""

def is_pvalue_column(col_name: str) -> bool:
    """Determine if a column is p-value (lower=better) or score (higher=better)."""

def load_truth_set(truth_path: str) -> tuple:
    """Load and validate truth set file. Returns (truth_df, truth_pos_df)."""

def compute_f1(precision: float, recall: float) -> float:
    """Compute F1 score from precision and recall."""

def bin_coverage(coverage_values: np.ndarray, bins: list = None) -> np.ndarray:
    """Bin coverage values into discrete categories."""

def generate_thresholds(scores: pd.Series, n_thresholds: int = 100) -> list:
    """Generate percentile-based threshold values."""
```

### 10.3 Script Functions (accuracy_benchmark.py)

```python
def detect_score_column(df: pd.DataFrame, tool_name: str = None) -> str:
    """Detect the most likely score column for ranking predictions."""

def compute_mcc(tp: int, fp: int, fn: int, tn: int) -> float:
    """Compute Matthews Correlation Coefficient."""

def compute_ranking_metrics(pred_df, truth_pos, truth_neg, score_col, window) -> tuple:
    """Compute AUPRC and AUROC given predictions with scores."""
```

### 10.4 Statistical Functions (benchmark_statistics.py)

```python
def bootstrap_metrics(metrics: np.ndarray, n_bootstrap: int = 1000,
                      alpha: float = 0.05) -> dict:
    """Compute bootstrap confidence intervals for metric arrays."""

def pairwise_wilcoxon(data: dict, metric: str) -> pd.DataFrame:
    """Perform pairwise Wilcoxon rank-sum tests between tools."""

def fdr_correction(p_values: np.ndarray, method: str = 'benjamini-hochberg') -> np.ndarray:
    """Apply FDR correction to p-values."""
```

### 10.5 Sensitivity Functions (benchmark_sensitivity.py)

```python
def stratify_by_coverage(df: pd.DataFrame, coverage_col: str,
                         bins: list) -> pd.DataFrame:
    """Add coverage bin column based on specified bins."""

def cross_validate_threshold(df: pd.DataFrame, score_col: str,
                             n_splits: int = 5) -> dict:
    """Cross-validation for threshold robustness analysis."""
```

---

## 11. Error Handling

### 11.1 Missing Truth Set
- If `benchmark.truth_set` is empty string → Skip accuracy benchmarking
- If file does not exist → Raise error in `accuracy_benchmark` rule

### 11.2 Missing Columns
- Truth set missing required columns → Python assertion error with clear message
- Tool output missing transcript/position → Warning logged, file skipped

### 11.3 Empty Results
- No predictions from any tool → Empty output file with header only
- No truth sites for modification type → Skipped for that type

### 11.4 Inferred Negatives
- When `label` column is missing, negatives are inferred from predictions
- `specificity` and `mcc` are set to `NaN` (not meaningful for inferred negatives)

### 11.5 Bootstrap Failures
- If bootstrap sampling fails (e.g., all NaN) → Return NaN for CIs with warning
- If too few samples for Wilcoxon → Report as "insufficient data"

---

## 12. Testing Strategy

### 12.1 Unit Tests
- Test `normalize_columns()` with various column name variants
- Test `detect_score_column()` with tool-specific dataframes
- Test `compute_mcc()` with edge cases (zeros, perfect scores)
- Test `is_pvalue_column()` with common column naming patterns

### 12.2 Integration Tests
- Run full benchmarking with sample truth set in `.test/`
- Verify output file format and metric calculations
- Test with multiple window values
- Verify statistical tests produce expected results

### 12.3 Test Data
- Synthetic truth set with known TP/FP/FN/TN counts
- Tool outputs with known score distributions
- Edge cases: empty files, single tool, single modification type
