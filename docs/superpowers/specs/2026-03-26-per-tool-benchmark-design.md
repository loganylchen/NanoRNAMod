# Per-Tool Benchmark System Design

**Date**: 2026-03-26
**Status**: Draft
**Author**: Claude (with user collaboration)

## 1. Overview

### 1.1 Problem Statement

The current benchmark system uses a single aggregated rule (`accuracy_benchmark`) that processes all modification detection tools together. This design has limitations:

- **No parallelization**: All tools must be processed sequentially in one rule
- **No incremental updates**: Adding a new tool requires re-running the entire benchmark
- **Poor modularity**: Tool-specific logic is mixed with aggregation logic

### 1.2 Proposed Solution

Refactor into per-tool benchmark rules that can run independently and in parallel, with a coverage-first approach for fair comparison across tools.

### 1.3 Design Goals

- **Parallelization**: Each tool's benchmark runs independently
- **Incremental updates**: Re-run only affected tools when data changes
- **Modularity**: Clear separation between tool-specific and aggregation logic
- **Fair comparison**: Coverage union ensures all tools evaluated on same sites
- **Native evaluation**: Tools also evaluated on only their covered sites

## 2. Coverage-First Design

### 2.1 Core Concept

For each comparison (native_sample vs control_sample):

1. **Compute coverage**: Each tool identifies which truth sites it can score
2. **Union coverage**: Take the union of all tools' covered sites
3. **Native evaluation**: Each tool evaluated only on sites it covers
4. **Fair evaluation**: All tools evaluated on union sites; missing predictions treated as negative (score = 0)

### 2.2 Coverage Definition

A site is "covered" by a tool if:
- The tool produces any prediction (with score) for that transcript/position
- No distinction between modification types - all treated as "modified base"

### 2.3 Evaluation Modes

| Mode | Evaluation Sites | Missing Predictions |
|------|-----------------|---------------------|
| Native | Tool's covered sites only | Excluded from calculation |
| Fair | Union of all tools' covered sites | Treated as negative (score = 0) |

## 3. Rule Structure

### 3.1 Coverage Rules

**Rule**: `benchmark_{tool}_coverage`

```
Input:  {project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv
        config["benchmark"]["truth_set"]

Output: {project}/results/benchmarks/coverage/{comparison}/{tool}_covered.tsv

Params: window (positional tolerance)
```

**Output format** (`{tool}_covered.tsv`):
```
transcript    position
ENST0001      1234
ENST0001      5678
...
```

**Logic**:
1. Load truth set
2. Load tool predictions
3. For each truth site, check if tool has prediction within window
4. Output covered sites

### 3.2 Coverage Union Rule

**Rule**: `benchmark_coverage_union`

```
Input:  {project}/results/benchmarks/coverage/{comparison}/{tool}_covered.tsv (for all tools)

Output: {project}/results/benchmarks/coverage/{comparison}/union.tsv

Params: None (aggregates all activated tools)
```

**Output format** (`union.tsv`):
```
transcript    position    covered_by_tools
ENST0001      1234        xpore,nanocompore,baleen
ENST0001      5678        xpore,baleen
...
```

**Logic**:
1. Collect all covered sites from each tool
2. Take union across all tools
3. Track which tools cover each site

### 3.3 Native Benchmark Rules

**Rule**: `benchmark_{tool}_native`

```
Input:  {project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv
        config["benchmark"]["truth_set"]

Output: {project}/results/benchmarks/native/{tool}/{comparison}/score_comparison.tsv
        {project}/results/benchmarks/native/{tool}/{comparison}/best_metrics.tsv
        {project}/results/benchmarks/native/{tool}/{comparison}/best_score.txt

Params: window (positional tolerance)
```

**Can run in parallel with coverage rules** - does not depend on union.

**Output format** (`score_comparison.tsv`):
```
score_column    is_pvalue    transform    auroc    prauc    best_threshold    best_threshold_original    f1    precision    recall
p_value         True         -log10       0.92     0.88     2.5               0.00316                    0.85  0.88         0.82
diff_mod        False        none         0.78     0.71     0.45              0.45                       0.71  0.75         0.68
...
```

**Output format** (`best_metrics.tsv`):
```
score_column    p_value_transform    auroc    prauc    best_threshold    best_threshold_original    f1    precision    recall
p_value         True                 0.92     0.88     2.5               0.00316                    0.85  0.88         0.82
```

**Output format** (`best_score.txt`):
```
p_value
```

### 3.4 Fair Benchmark Rules

**Rule**: `benchmark_{tool}_fair`

```
Input:  {project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv
        {project}/results/benchmarks/coverage/{comparison}/union.tsv
        config["benchmark"]["truth_set"]

Output: {project}/results/benchmarks/fair/{tool}/{comparison}/score_comparison.tsv
        {project}/results/benchmarks/fair/{tool}/{comparison}/best_metrics.tsv
        {project}/results/benchmarks/fair/{tool}/{comparison}/best_score.txt

Params: window (positional tolerance)
```

**Depends on**: `benchmark_coverage_union` (must wait for all coverage rules)

**Output format**: Same as native benchmark rules

**Logic**:
1. Load union coverage sites
2. Load tool predictions
3. For sites in union but not covered by tool: assign score = 0 (negative prediction)
4. Evaluate all candidate score columns
5. Report best column by AUROC

## 4. Score Column Handling

### 4.1 Score Column Detection

Each tool has multiple columns that could serve as scores. The system automatically detects and evaluates all candidates.

**Per-tool score column priority** (from existing `detect_score_column()`):
```python
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
```

### 4.2 Score Column Processing

For each candidate score column:

```python
def process_score_column(df, score_col):
    # 1. Type conversion
    scores = pd.to_numeric(df[score_col], errors='coerce')

    # 2. Handle special values
    scores = scores.replace([np.inf], 1e10)
    scores = scores.replace([-np.inf], -1e10)
    # NA/NaN and strings already coerced to NaN

    # 3. P-value detection and transformation
    is_pvalue = is_pvalue_column(score_col)
    if is_pvalue:
        scores = -np.log10(scores)

    return scores, is_pvalue
```

**P-value detection**:
```python
def is_pvalue_column(col_name):
    pvalue_keywords = ['p_value', 'pvalue', 'pval', 'p.value']
    col_lower = col_name.lower().replace(' ', '_')
    return any(kw in col_lower for kw in pvalue_keywords)
```

### 4.3 Threshold Optimization

For each score column, find the optimal threshold that maximizes F1:

```python
def find_optimal_threshold(scores, labels, is_pvalue):
    """
    Args:
        scores: Array of scores (higher = better after transformation)
        labels: Binary labels (1 = modified, 0 = unmodified)
        is_pvalue: Whether this was -log10 transformed

    Returns:
        dict with threshold, metrics
    """
    # Generate thresholds from score distribution
    thresholds = np.linspace(scores.min(), scores.max(), 100)

    best_f1 = 0
    best_threshold = None
    best_metrics = None

    for thresh in thresholds:
        predictions = (scores >= thresh).astype(int)
        metrics = compute_metrics(predictions, labels)

        if metrics['f1'] > best_f1:
            best_f1 = metrics['f1']
            best_threshold = thresh
            best_metrics = metrics

    # Convert threshold back to original scale
    if is_pvalue:
        original_threshold = 10 ** (-best_threshold)
    else:
        original_threshold = best_threshold

    return {
        'threshold': best_threshold,
        'threshold_original': original_threshold,
        'f1': best_metrics['f1'],
        'precision': best_metrics['precision'],
        'recall': best_metrics['recall'],
    }
```

### 4.4 Metric Computation

For each score column:

```python
def compute_all_metrics(scores, labels):
    """
    Compute all metrics for a score column.

    Args:
        scores: Array of scores (higher = better)
        labels: Binary labels

    Returns:
        dict with auroc, prauc, f1, precision, recall, best_threshold, etc.
    """
    from sklearn.metrics import roc_auc_score, average_precision_score

    auroc = roc_auc_score(labels, scores)
    prauc = average_precision_score(labels, scores)

    # Find optimal threshold
    threshold_result = find_optimal_threshold(scores, labels, is_pvalue=False)

    return {
        'auroc': auroc,
        'prauc': prauc,
        'f1': threshold_result['f1'],
        'precision': threshold_result['precision'],
        'recall': threshold_result['recall'],
        'best_threshold': threshold_result['threshold'],
        'best_threshold_original': threshold_result['threshold_original'],
    }
```

### 4.5 Score Comparison Output

All score columns are evaluated and ranked by AUROC. The best score column is identified, but all results are retained.

**`score_comparison.tsv`** columns:
| Column | Description |
|--------|-------------|
| score_column | Name of the score column |
| is_pvalue | Whether this was -log10 transformed |
| transform | Transformation applied ('-log10' or 'none') |
| auroc | Area under ROC curve |
| prauc | Area under PR curve |
| best_threshold | Threshold in transformed space |
| best_threshold_original | Threshold in original scale |
| f1 | F1 score at optimal threshold |
| precision | Precision at optimal threshold |
| recall | Recall at optimal threshold |

## 5. Aggregation Rules

### 5.1 Coverage Union Aggregation

**Rule**: `aggregate_coverage_union`

Already covered by `benchmark_coverage_union` - produces per-comparison union files.

### 5.2 Tool-Level Aggregation

**Rule**: `aggregate_tool_results`

```
Input:  {project}/results/benchmarks/native/{tool}/{comparison}/*.tsv
        {project}/results/benchmarks/fair/{tool}/{comparison}/*.tsv

Output: {project}/results/benchmarks/aggregated/{tool}/native_summary.tsv
        {project}/results/benchmarks/aggregated/{tool}/fair_summary.tsv
```

**Output format** (`native_summary.tsv`):
```
tool    best_score    mean_auroc    std_auroc    mean_prauc    std_prauc    mean_f1    std_f1    n_comparisons
xpore   p_value       0.89          0.05         0.85          0.06         0.82       0.04      3
```

### 5.3 Cross-Tool Comparison

**Rule**: `aggregate_benchmark_comparison`

```
Input:  {project}/results/benchmarks/aggregated/{tool}/*_summary.tsv (for all tools)

Output: {project}/results/benchmarks/summary/by_comparison.tsv
        {project}/results/benchmarks/summary/by_tool.tsv
        {project}/results/benchmarks/summary/best_scores.tsv
```

**Output format** (`by_comparison.tsv`):
```
comparison    xpore_auroc    nanocompore_auroc    baleen_auroc    ...    best_tool
A_F           0.92           0.88                 0.85            ...    xpore
B_G           0.85           0.90                 0.82            ...    nanocompore
```

**Output format** (`by_tool.tsv`):
```
tool         mean_auroc    std_auroc    mean_prauc    mean_f1    n_comparisons
xpore        0.88          0.05         0.84          0.80       3
nanocompore  0.87          0.04         0.83          0.79       3
baleen       0.82          0.06         0.78          0.75       3
```

**Output format** (`best_scores.tsv`):
```
tool         best_score    n_comparisons_best
xpore        p_value       2
nanocompore  logit_pvalue  1
```

## 6. Directory Structure

```
{project}/results/benchmarks/
├── coverage/
│   └── {comparison}/
│       ├── {tool}_covered.tsv
│       └── union.tsv
├── native/
│   └── {tool}/
│       └── {comparison}/
│           ├── score_comparison.tsv
│           ├── best_metrics.tsv
│           └── best_score.txt
├── fair/
│   └── {tool}/
│       └── {comparison}/
│           ├── score_comparison.tsv
│           ├── best_metrics.tsv
│           └── best_score.txt
├── aggregated/
│   └── {tool}/
│       ├── native_summary.tsv
│       └── fair_summary.tsv
└── summary/
    ├── by_comparison.tsv
    ├── by_tool.tsv
    └── best_scores.tsv
```

## 7. Execution DAG

```
                    ┌─────────────────────────────────────┐
                    │  benchmark_{tool}_coverage          │
                    │  (parallel per tool)                │
                    └──────────────┬──────────────────────┘
                                   │
                    ┌──────────────▼──────────────────────┐
                    │  benchmark_coverage_union           │
                    │  (waits for all coverage rules)     │
                    └──────────────┬──────────────────────┘
                                   │
          ┌────────────────────────┼────────────────────────┐
          │                        │                        │
          ▼                        ▼                        ▼
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│ benchmark_xpore │    │benchmark_nanocom│    │ benchmark_baleen│
│    _native      │    │   pore_native   │    │    _native      │
└─────────────────┘    └─────────────────┘    └─────────────────┘
          │                        │                        │
          ▼                        ▼                        ▼
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│ benchmark_xpore │    │benchmark_nanocom│    │ benchmark_baleen│
│    _fair        │    │   pore_fair     │    │    _fair        │
└─────────────────┘    └─────────────────┘    └─────────────────┘
          │                        │                        │
          └────────────────────────┼────────────────────────┘
                                   │
                    ┌──────────────▼──────────────────────┐
                    │  aggregate_tool_results             │
                    │  (per tool, parallel)               │
                    └──────────────┬──────────────────────┘
                                   │
                    ┌──────────────▼──────────────────────┐
                    │  aggregate_benchmark_comparison     │
                    │  (final summary)                    │
                    └─────────────────────────────────────┘
```

**Parallelization opportunities**:
- All `benchmark_{tool}_coverage` rules run in parallel
- All `benchmark_{tool}_native` rules run in parallel
- All `benchmark_{tool}_fair` rules run in parallel
- All `aggregate_tool_results` rules run in parallel

## 8. Implementation Notes

### 8.1 Backward Compatibility

The new structure should coexist with the existing benchmark system during transition. Consider:

1. Keep existing `accuracy_summary.tsv` outputs for downstream rules
2. Add new outputs in parallel
3. Gradually migrate downstream rules to use new structure

### 8.2 Configuration

Add to `config/config.yaml`:

```yaml
benchmark:
  truth_set: "path/to/truth.tsv"
  window: 0  # positional tolerance
  criterion: "f1"  # threshold optimization criterion
  modes:  # which modes to run
    - native
    - fair
```

### 8.3 Rule Organization

Create new rule files:
- `workflow/rules/benchmark_coverage.smk` - coverage and union rules
- `workflow/rules/benchmark_per_tool.smk` - native and fair per-tool rules
- `workflow/rules/benchmark_aggregation.smk` - aggregation rules

Or keep in existing `benchmark_accuracy.smk` with clear sections.

### 8.4 Script Organization

Create new scripts:
- `workflow/scripts/benchmark_coverage.py` - coverage computation
- `workflow/scripts/benchmark_per_tool.py` - per-tool native/fair evaluation
- `workflow/scripts/benchmark_aggregation.py` - result aggregation

## 9. Testing Strategy

### 9.1 Unit Tests

- Score column detection and transformation
- Threshold optimization
- Coverage union computation
- Metric computation

### 9.2 Integration Tests

- Run full pipeline on test dataset
- Verify output format compliance
- Check parallelization works correctly

### 9.3 Validation

- Compare results with existing benchmark system
- Verify AUROC/PRAUC calculations match sklearn
- Confirm threshold conversion is correct

## 10. Future Enhancements

### 10.1 Additional Metrics

- MCC (Matthews Correlation Coefficient)
- Balanced accuracy
- Specificity
- NPV (Negative Predictive Value)

### 10.2 Visualization

- Per-tool ROC/PR curves
- Coverage heatmaps
- Score distribution plots
- Threshold sensitivity analysis

### 10.3 Statistical Analysis

- Bootstrap confidence intervals
- Significance testing between tools
- Effect size calculations
