# Integrated Benchmark System Design

**Date**: 2026-03-25
**Status**: Draft for Review
**Goal**: Integrate comprehensive, publication-ready benchmarking into the NanoRNAMod Snakemake workflow

---

## 1. Overview

### 1.1 Problem Statement

The current benchmark system provides solid foundational metrics but lacks:
- Statistical rigor (no confidence intervals, no significance testing)
- Publication-ready figure aesthetics (Nature/Cell/Science standards)
- Systematic sensitivity analysis
- Automated figure legend and caption generation
- Supplementary materials organization

### 1.2 Design Goals

1. **Workflow Integration**: Every dataset run through NanoRNAMod automatically produces complete benchmark results
2. **Statistical Rigor**: Bootstrap CIs, paired significance tests, FDR correction
3. **Publication Quality**: Nature/Cell/Science figure standards with proper legends
4. **Reproducibility**: All plotting data preserved, all parameters logged
5. **Extensibility**: Easy to add new tools or metrics

### 1.3 Scope

**In Scope**:
- 9 comparison-based tools: xpore, nanocompore, baleen, differr, drummer, eligos2, epinano, psipore, pybaleen
- Multi-threshold evaluation (50 thresholds default)
- Score column optimization per tool
- Two negative control strategies (k-mer, same-base)
- Ground truth: Literature-curated sites (m6A, m5C, ψ, m1A, etc.)

**Out of Scope**:
- Per-sample tools (tandemmod, directrm, m6atm, rnano) - future phase
- Real-time benchmarking dashboards
- Cloud-based result aggregation

---

## 2. Architecture

### 2.1 Directory Structure

```
{project}/results/benchmarks/
├── .benchmark_complete              # Touch file for downstream rules
├── accuracy_summary.tsv             # Aggregated metrics (existing)
├── accuracy_summary_overall.tsv     # Overall metrics (existing)
├── accuracy_summary_by_comparison.tsv
├── accuracy_summary_by_negative_type.tsv
├── accuracy_summary_kmer_negatives.tsv
├── accuracy_summary_same_base_negatives.tsv
├── kmer_negative_sites.tsv
├── same_base_negative_sites.tsv
│
├── statistics/                      # NEW: Statistical analysis
│   ├── bootstrap_ci.tsv             # Bootstrap confidence intervals
│   ├── significance_tests.tsv       # Paired Wilcoxon + permutation tests
│   ├── fdr_corrected.tsv            # FDR-corrected p-values
│   └── effect_sizes.tsv             # Cohen's d, odds ratios
│
├── sensitivity/                     # NEW: Sensitivity analysis
│   ├── coverage_analysis.tsv        # Depth vs accuracy
│   ├── score_distribution.tsv       # Score distribution by truth status
│   └── threshold_robustness.tsv     # Threshold stability analysis
│
├── figures/                         # NEW: Publication-ready figures
│   ├── main/                        # Main text figures
│   │   ├── fig1_pr_curves.pdf
│   │   ├── fig2_roc_curves.pdf
│   │   ├── fig3_heatmap.pdf
│   │   ├── fig4_resource_comparison.pdf
│   │   └── fig5_optimal_thresholds.pdf
│   ├── supplementary/               # Supplementary figures
│   │   ├── sfig1_per_comparison.pdf
│   │   ├── sfig2_score_distributions.pdf
│   │   └── ...
│   └── legends/                     # Auto-generated figure legends
│       ├── fig1_legend.md
│       ├── fig2_legend.md
│       └── ...
│
├── data/                            # NEW: Preserved plotting data
│   ├── pr_curve_data.tsv
│   ├── roc_curve_data.tsv
│   ├── heatmap_data.tsv
│   └── ...
│
├── supplementary/                   # NEW: Supplementary materials
│   ├── supplementary_methods.md
│   ├── supplementary_tables.tsv
│   └── supplementary_data.tsv
│
├── viz/                             # Existing: HTML report
│   └── benchmark_report.html
│
├── threshold_evaluation.tsv         # Existing
├── optimal_thresholds_detailed.tsv  # Existing
├──── optimal_score_per_tool.tsv     # Existing
├── score_distributions.tsv          # Existing
├── detailed_predictions.tsv         # Existing
├── resource_summary.tsv             # Existing
├── resource_by_tool.tsv             # Existing
└── benchmark_report.pdf             # Existing
```

### 2.2 Rule Dependency Graph

```
Tool Results (modetect_*)
        │
        ▼
accuracy_benchmark ─────────────────────────────────────────┐
        │                                                    │
        ├─► benchmark_kmer_negatives                         │
        ├─► benchmark_same_base_negatives                    │
        ├─► benchmark_multithreshold                         │
        └─► benchmark_score_optimization                     │
        │                                                    │
        ▼                                                    │
benchmark_statistics (NEW) ◄─────────────────────────────────┤
        │                                                    │
        ├─► Bootstrap CI computation                         │
        ├─► Paired significance tests                        │
        └─► FDR correction                                   │
        │                                                    │
        ▼                                                    │
benchmark_sensitivity (NEW)                                  │
        │                                                    │
        ├─► Coverage depth analysis                          │
        ├─► Score distribution analysis                      │
        └─► Threshold robustness analysis                    │
        │                                                    │
        ▼                                                    │
benchmark_visualization (existing, enhanced)                  │
        │                                                    │
        ▼                                                    │
benchmark_figures (NEW)                                      │
        │                                                    │
        ├─► Main figures with Nature theme                   │
        ├─► Supplementary figures                            │
        └─► Auto-generated legends                           │
        │                                                    │
        ▼                                                    │
benchmark_supplementary (NEW)                                │
        │                                                    │
        ├─► Methods documentation                            │
        ├─► Supplementary tables                             │
        └─► Data exports                                     │
        │                                                    │
        ▼                                                    │
benchmark_report (enhanced PDF)                              │
```

---

## 3. Component Design

### 3.1 Statistical Analysis Module

**File**: `workflow/scripts/benchmark_statistics.py`

**Functions**:

```python
def compute_bootstrap_ci(df, metric_col, n_bootstrap=1000, ci=95):
    """Compute bootstrap confidence intervals for a metric."""
    # Stratified bootstrap by comparison
    # Returns: (lower, upper, point_estimate)

def paired_wilcoxon_test(df, tool1, tool2, metric_col):
    """Paired Wilcoxon signed-rank test between tools."""
    # Returns: statistic, p-value, effect_size

def permutation_test(df, tool1, tool2, metric_col, n_perm=10000):
    """Permutation test for robust significance."""
    # Returns: p-value, null_distribution

def fdr_correction(p_values, method='benjamini-hochberg'):
    """Apply FDR correction to p-values."""
    # Returns: adjusted_p_values

def compute_effect_sizes(df, tool1, tool2, metric_col):
    """Compute Cohen's d and odds ratio."""
    # Returns: cohens_d, odds_ratio, ci
```

**Output Format** (`statistics/bootstrap_ci.tsv`):
```
tool	modification_type	metric	point_estimate	ci_lower	ci_upper	n_bootstrap
xpore	m6A	f1	0.72	0.68	0.76	1000
xpore	m6A	precision	0.81	0.77	0.85	1000
...
```

**Output Format** (`statistics/significance_tests.tsv`):
```
tool1	tool2	modification_type	metric	test_type	statistic	p_value	adj_p_value	effect_size
xpore	baleen	m6A	f1	wilcoxon	1234.5	0.023	0.045	0.34
xpore	nanocompore	m6A	f1	permutation	-	0.031	0.058	0.28
...
```

### 3.2 Sensitivity Analysis Module

**File**: `workflow/scripts/benchmark_sensitivity.py`

**Functions**:

```python
def analyze_coverage_effect(df, coverage_bins=None):
    """Analyze accuracy vs coverage depth."""
    # Default bins: [0, 10, 20, 50, 100, 200, 500, inf]
    # Returns: accuracy by coverage bin

def analyze_score_distribution(df, tool, score_col):
    """Analyze score distribution by truth status."""
    # Returns: distribution statistics, KS test

def analyze_threshold_robustness(df, n_splits=5):
    """Assess threshold stability via cross-validation."""
    # Returns: threshold variance across splits
```

**Output Format** (`sensitivity/coverage_analysis.tsv`):
```
tool	coverage_bin	n_sites	precision	recall	f1
xpore	[0,10)	150	0.45	0.32	0.37
xpore	[10,20)	230	0.58	0.51	0.54
...
```

### 3.3 Publication Visualization Module

**File**: `workflow/scripts/benchmark_figures.py`

**Theme Configuration**:
```python
NATURE_THEME = {
    'font_family': 'Helvetica',
    'font_size': 7,  # Nature standard
    'figsize': (3.5, 2.5),  # Single column width in inches at 300dpi
    'dpi': 300,
    'linewidth': 0.5,
    'markersize': 3,
    'colors': {
        'primary': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
        'positive': '#2ca02c',
        'negative': '#d62728',
        'neutral': '#7f7f7f'
    },
    'grid': {'visible': True, 'alpha': 0.3, 'linewidth': 0.5},
    'spines': {'visible': ['left', 'bottom'], 'linewidth': 0.5},
}
```

**Figure Functions**:

```python
def fig1_pr_curves(df, output_path):
    """Main figure 1: PR curves with CIs for all tools."""
    # - Multi-panel layout (2x5 for 9 tools)
    # - Shaded CI regions
    # - AUROC/AUPRC in legend
    # - Nature theme styling

def fig2_roc_curves(df, output_path):
    """Main figure 2: ROC curves with CIs for all tools."""
    # Same layout as PR curves

def fig3_performance_heatmap(df, output_path):
    """Main figure 3: Heatmap of metrics by tool and modification type."""
    # - F1 score as primary metric
    # - Hierarchical clustering
    # - Colorblind-friendly palette

def fig4_resource_comparison(df, output_path):
    """Main figure 4: Runtime and memory comparison."""
    # - Bar plot with error bars
    # - Log scale for memory
    # - Tool names on x-axis

def fig5_optimal_thresholds(df, output_path):
    """Main figure 5: Optimal threshold distribution by tool."""
    # - Violin or box plots
    # - Individual points
    # - Statistical annotations
```

### 3.4 Legend Generation Module

**File**: `workflow/scripts/benchmark_legends.py`

**Template**:
```python
LEGEND_TEMPLATES = {
    'fig1_pr_curves': """
**Figure 1: Precision-Recall Analysis of RNA Modification Detection Tools.**

Precision-recall curves for {n_tools} comparison-based tools evaluated on {n_sites} literature-curated modification sites across {n_mods} modification types (m6A, m5C, ψ, m1A). Shaded regions indicate 95% bootstrap confidence intervals (n=1,000 iterations). Area under the PR curve (AUPRC) values are shown in parentheses. Tools are ordered by overall F1 score. Dashed diagonal line indicates random classifier baseline.
""",
    # ... more templates
}

def generate_legend(figure_name, context):
    """Generate figure legend with populated context."""
    template = LEGEND_TEMPLATES[figure_name]
    return template.format(**context)
```

### 3.5 Supplementary Materials Module

**File**: `workflow/scripts/benchmark_supplementary.py`

**Contents**:

1. **Supplementary Methods** (`supplementary/supplementary_methods.md`):
   - Benchmark dataset description
   - Tool version and parameters
   - Statistical methods
   - Figure generation code

2. **Supplementary Tables** (`supplementary/supplementary_tables.tsv`):
   - Complete metric tables
   - Per-comparison breakdown
   - Statistical test results

3. **Supplementary Data** (`supplementary/supplementary_data.tsv`):
   - All site-level predictions
   - Truth annotations
   - Coverage information

---

## 4. Snakemake Rule Design

### 4.1 New Rule: benchmark_statistics

```python
rule benchmark_statistics:
    """Compute bootstrap CIs and significance tests."""
    input:
        summary="{project}/results/benchmarks/accuracy_summary.tsv",
        by_comparison="{project}/results/benchmarks/accuracy_summary_by_comparison.tsv",
    output:
        ci="{project}/results/benchmarks/statistics/bootstrap_ci.tsv",
        significance="{project}/results/benchmarks/statistics/significance_tests.tsv",
        fdr="{project}/results/benchmarks/statistics/fdr_corrected.tsv",
        effects="{project}/results/benchmarks/statistics/effect_sizes.tsv",
    params:
        n_bootstrap=config.get("benchmark", {}).get("n_bootstrap", 1000),
        alpha=config.get("benchmark", {}).get("alpha", 0.05),
    resources:
        mem_mb=1024 * 8,
    threads: 4
    priority: 18
    log:
        "logs/{project}/benchmark_statistics/statistics.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_statistics.py"
```

### 4.2 New Rule: benchmark_sensitivity

```python
rule benchmark_sensitivity:
    """Sensitivity analysis for coverage and threshold robustness."""
    input:
        predictions="{project}/results/benchmarks/detailed_predictions.tsv",
        truth_set=config["benchmark"]["truth_set"],
    output:
        coverage="{project}/results/benchmarks/sensitivity/coverage_analysis.tsv",
        score_dist="{project}/results/benchmarks/sensitivity/score_distribution.tsv",
        threshold_robust="{project}/results/benchmarks/sensitivity/threshold_robustness.tsv",
    params:
        coverage_bins=config.get("benchmark", {}).get("coverage_bins", [0, 10, 20, 50, 100, 200, 500]),
        n_splits=config.get("benchmark", {}).get("n_splits", 5),
    resources:
        mem_mb=1024 * 8,
    threads: 2
    priority: 17
    log:
        "logs/{project}/benchmark_sensitivity/sensitivity.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_sensitivity.py"
```

### 4.3 New Rule: benchmark_figures

```python
rule benchmark_figures:
    """Generate publication-ready figures with Nature theme."""
    input:
        summary="{project}/results/benchmarks/accuracy_summary.tsv",
        ci="{project}/results/benchmarks/statistics/bootstrap_ci.tsv",
        significance="{project}/results/benchmarks/statistics/significance_tests.tsv",
        coverage="{project}/results/benchmarks/sensitivity/coverage_analysis.tsv",
        resources="{project}/results/benchmarks/resource_summary.tsv",
        thresholds="{project}/results/benchmarks/threshold_evaluation.tsv",
    output:
        main=expand("{{project}}/results/benchmarks/figures/main/fig{n}.pdf", n=[1,2,3,4,5]),
        supplementary=directory("{project}/results/benchmarks/figures/supplementary/"),
        legends=expand("{{project}}/results/benchmarks/figures/legends/fig{n}_legend.md", n=[1,2,3,4,5]),
        data=directory("{project}/results/benchmarks/data/"),
    params:
        theme=config.get("benchmark", {}).get("figure_theme", "nature"),
        dpi=config.get("benchmark", {}).get("dpi", 300),
    resources:
        mem_mb=1024 * 4,
    threads: 1
    priority: 40
    log:
        "logs/{project}/benchmark_figures/figures.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_figures.py"
```

### 4.4 New Rule: benchmark_supplementary

```python
rule benchmark_supplementary:
    """Generate supplementary materials."""
    input:
        summary="{project}/results/benchmarks/accuracy_summary.tsv",
        predictions="{project}/results/benchmarks/detailed_predictions.tsv",
        statistics=expand("{{project}}/results/benchmarks/statistics/{f}",
                         f=["bootstrap_ci.tsv", "significance_tests.tsv", "effect_sizes.tsv"]),
    output:
        methods="{project}/results/benchmarks/supplementary/supplementary_methods.md",
        tables="{project}/results/benchmarks/supplementary/supplementary_tables.tsv",
        data="{project}/results/benchmarks/supplementary/supplementary_data.tsv",
    params:
        tools=lambda: [t for t in config["tools"] if config["tools"][t]["activate"]],
        truth_set=config["benchmark"]["truth_set"],
    resources:
        mem_mb=1024 * 2,
    threads: 1
    priority: 45
    log:
        "logs/{project}/benchmark_supplementary/supplementary.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_supplementary.py"
```

### 4.5 Modified Rule: get_final_output

Add new outputs to `workflow/rules/common.smk`:

```python
def get_final_output():
    # ... existing code ...

    if config.get("benchmark", {}).get("truth_set", ""):
        # ... existing benchmark outputs ...

        # NEW: Statistical analysis
        final_output += [f"{RESULT_ROOT}/benchmarks/statistics/bootstrap_ci.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/statistics/significance_tests.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/statistics/fdr_corrected.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/statistics/effect_sizes.tsv"]

        # NEW: Sensitivity analysis
        final_output += [f"{RESULT_ROOT}/benchmarks/sensitivity/coverage_analysis.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/sensitivity/score_distribution.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/sensitivity/threshold_robustness.tsv"]

        # NEW: Publication figures
        final_output += expand(f"{RESULT_ROOT}/benchmarks/figures/main/fig{{n}}.pdf", n=[1,2,3,4,5])
        final_output += [f"{RESULT_ROOT}/benchmarks/figures/supplementary/"]
        final_output += expand(f"{RESULT_ROOT}/benchmarks/figures/legends/fig{{n}}_legend.md", n=[1,2,3,4,5])
        final_output += [f"{RESULT_ROOT}/benchmarks/data/"]

        # NEW: Supplementary materials
        final_output += [f"{RESULT_ROOT}/benchmarks/supplementary/supplementary_methods.md"]
        final_output += [f"{RESULT_ROOT}/benchmarks/supplementary/supplementary_tables.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/supplementary/supplementary_data.tsv"]

    return final_output
```

---

## 5. Configuration Schema Updates

### 5.1 config/config.yaml additions

```yaml
benchmark:
  truth_set: "config/benchmark_ecoli_rRNA.tsv"
  window: [0, 1, 2, 5]
  n_thresholds: 50
  custom_thresholds: []

  # NEW: Statistical analysis
  n_bootstrap: 1000
  alpha: 0.05
  fdr_method: "benjamini-hochberg"

  # NEW: Sensitivity analysis
  coverage_bins: [0, 10, 20, 50, 100, 200, 500]
  n_splits: 5

  # NEW: Figure settings
  figure_theme: "nature"  # Options: nature, cell, science
  dpi: 300
  fig_format: "pdf"  # Options: pdf, png, svg

  # NEW: Supplementary materials
  include_supplementary: true
  data_export_format: "tsv"  # Options: tsv, csv, parquet
```

### 5.2 Schema validation updates

Add to `schemas/config.schema.yaml`:

```yaml
benchmark:
  type: object
  properties:
    truth_set:
      type: string
    window:
      type: array
      items:
        type: integer
    n_thresholds:
      type: integer
      minimum: 10
    custom_thresholds:
      type: array
    # NEW properties
    n_bootstrap:
      type: integer
      minimum: 100
      default: 1000
    alpha:
      type: number
      minimum: 0
      maximum: 1
      default: 0.05
    fdr_method:
      type: string
      enum: ["benjamini-hochberg", "bonferroni", "holm"]
      default: "benjamini-hochberg"
    coverage_bins:
      type: array
      items:
        type: integer
    n_splits:
      type: integer
      minimum: 2
      default: 5
    figure_theme:
      type: string
      enum: ["nature", "cell", "science"]
      default: "nature"
    dpi:
      type: integer
      minimum: 72
      default: 300
    fig_format:
      type: string
      enum: ["pdf", "png", "svg"]
      default: "pdf"
    include_supplementary:
      type: boolean
      default: true
    data_export_format:
      type: string
      enum: ["tsv", "csv", "parquet"]
      default: "tsv"
```

---

## 6. Implementation Plan

### Phase 1: Statistical Foundation (Priority: High)
1. Create `benchmark_statistics.py` with bootstrap CI and significance tests
2. Add `benchmark_statistics` rule to `benchmark_viz.smk`
3. Update `get_final_output()` with new outputs
4. Test with existing benchmark data

### Phase 2: Sensitivity Analysis (Priority: High)
1. Create `benchmark_sensitivity.py` with coverage and threshold analysis
2. Add `benchmark_sensitivity` rule
3. Integrate coverage data from alignment steps
4. Update final outputs

### Phase 3: Publication Figures (Priority: High)
1. Create `benchmark_figures.py` with Nature theme
2. Create `benchmark_legends.py` for auto-generated legends
3. Add `benchmark_figures` rule
4. Ensure all plotting data is preserved

### Phase 4: Supplementary Materials (Priority: Medium)
1. Create `benchmark_supplementary.py`
2. Add `benchmark_supplementary` rule
3. Generate methods documentation template
4. Create supplementary data exports

### Phase 5: Integration and Testing (Priority: High)
1. Update schema validation
2. Run full workflow test in `.test/`
3. Validate figure output meets Nature standards
4. Documentation updates

---

## 7. Dependencies

### 7.1 New Python packages

Add to `workflow/envs/python3.yaml`:

```yaml
dependencies:
  # ... existing ...
  - scipy>=1.10.0      # For statistical tests
  - statsmodels>=0.14  # For FDR correction
  - seaborn>=0.12      # For enhanced visualizations
  - matplotlib>=3.7    # For publication figures
```

### 7.2 Optional: R integration

For advanced statistical visualizations, consider adding R support:
- `ggplot2` for additional figure options
- `ggpubr` for publication-ready themes
- `rstatix` for statistical annotations

---

## 8. Testing Strategy

### 8.1 Unit tests

- Bootstrap CI: Known distribution (normal, binomial)
- Significance tests: Known differences
- FDR correction: Known p-value sets
- Figure generation: Visual inspection + size validation

### 8.2 Integration tests

- Full workflow run in `.test/`
- Verify all outputs generated
- Check figure dimensions and DPI
- Validate legend content

### 8.3 Validation criteria

- Bootstrap CIs contain true value 95% of time
- Significance tests detect known differences
- Figures meet Nature specifications (single column: 88mm, double column: 180mm)
- All data files loadable and complete

---

## 9. Future Extensions

1. **Per-sample tool benchmarking**: Extend to tandemmod, directrm, m6atm, rnano
2. **Cross-dataset comparison**: Aggregate benchmarks across multiple runs
3. **Interactive visualization**: Shiny/streamlit dashboard for exploration
4. **Automated report generation**: Word/LaTeX integration for manuscript prep
5. **Version tracking**: Track benchmark performance over code changes

---

## 10. Acceptance Criteria

- [ ] Bootstrap CIs computed for all metrics with configurable iterations
- [ ] Paired significance tests between all tool pairs with FDR correction
- [ ] Sensitivity analysis shows coverage/threshold effects
- [ ] All main figures generated in Nature format (PDF, 300 DPI)
- [ ] Auto-generated figure legends with populated statistics
- [ ] Supplementary materials organized and complete
- [ ] All plotting data preserved in TSV format
- [ ] Full workflow runs successfully in `.test/`
- [ ] Documentation updated with new benchmark features
