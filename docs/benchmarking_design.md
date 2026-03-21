# NanoRNAMod Benchmarking Module Design Document

## 1. Overview

### 1.1 Purpose
The benchmarking module provides standardized accuracy and performance evaluation for RNA modification detection tools integrated into NanoRNAMod. It enables:
- **Accuracy benchmarking**: Compare tool predictions against known ground truth modification sites
- **Resource benchmarking**: Track computational resources (CPU, memory, I/O) consumed by each tool
- **Cross-tool comparison**: Generate unified reports comparing all activated tools

### 1.2 Current Implementation Status

The following components are already implemented:
- `workflow/rules/benchmark_accuracy.smk` - Accuracy evaluation rule
- `workflow/rules/benchmark_report.smk` - Resource aggregation rule
- `workflow/scripts/accuracy_benchmark.py` - Core accuracy metrics computation
- `workflow/scripts/aggregate_benchmarks.py` - Resource benchmark aggregation
- Configuration schema for `benchmark.truth_set` and `benchmark.window`

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

---

## 4. Output File Structure

### 4.1 Directory Layout

```
{project}/results/benchmarks/
├── accuracy_summary.tsv           # Per-modification_type metrics
├── accuracy_summary_overall.tsv   # Aggregated across all mod types
└── resource_summary.tsv           # CPU/memory/IO benchmarks

logs/{project}/
├── accuracy_benchmark/
│   └── accuracy.log
└── aggregate_benchmarks/
    └── aggregate.log

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

### 4.4 resource_summary.tsv Columns

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
  window: [0, 1, 2, 5]  # Positional tolerance windows (list of integers)
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
| `workflow/rules/benchmark_report.smk` | `aggregate_benchmarks` | Aggregate resource benchmarks |

### 6.2 Existing Scripts

| Script File | Purpose |
|-------------|---------|
| `workflow/scripts/accuracy_benchmark.py` | Core accuracy computation logic |
| `workflow/scripts/aggregate_benchmarks.py` | Resource benchmark aggregation |

### 6.3 Conda Environments

| Environment File | Dependencies |
|------------------|--------------|
| `workflow/envs/pandas.yaml` | pandas, numpy, scikit-learn |

---

## 7. Proposed Enhancements

### 7.1 Per-Tool Detailed Reports

**New output**: `{project}/results/benchmarks/{tool}_detailed.tsv`

Contains site-by-site comparison:
```tsv
transcript	position	predicted	truth_match	distance_to_truth	score
```

**New script**: `workflow/scripts/detailed_benchmark.py`

### 7.2 Visualization Module

**New outputs**:
```
{project}/results/benchmarks/
├── figures/
│   ├── precision_recall_curve.svg
│   ├── roc_curve.svg
│   ├── f1_by_window.svg
│   └── resource_comparison.svg
└── benchmark_report.html
```

**New rule**: `benchmark_visualization` in `workflow/rules/benchmark_viz.smk`
**New script**: `workflow/scripts/benchmark_plots.py`

### 7.3 Threshold Optimization

**Feature**: Automatically find optimal score thresholds for each tool

**New output**: `{project}/results/benchmarks/optimal_thresholds.tsv`
```tsv
tool	modification_type	optimal_threshold	optimal_f1
```

### 7.4 Multi-Truth-Set Support

**Configuration**:
```yaml
benchmark:
  truth_sets:
    - name: "ecoli_rrna"
      path: "config/benchmark_ecoli_rRNA.tsv"
    - name: "hekar_m6a"
      path: "config/benchmark_HEKAR_m6A.tsv"
```

### 7.5 Stratified Analysis

**Feature**: Break down metrics by:
- Transcript
- Modification type
- Position context (e.g., DRACH motif for m6A)

---

## 8. Implementation Checklist

### Phase 1: Core Enhancements
- [ ] Add `modification_name` column support to truth set parsing
- [ ] Implement per-tool detailed reports
- [ ] Add threshold optimization logic
- [ ] Improve score column detection with tool-specific mappings

### Phase 2: Visualization
- [ ] Create `benchmark_viz.smk` rule file
- [ ] Implement `benchmark_plots.py` with matplotlib/seaborn
- [ ] Generate HTML benchmark report

### Phase 3: Advanced Features
- [ ] Multi-truth-set support
- [ ] Stratified analysis by transcript/motif
- [ ] Confidence interval calculation (bootstrapping)

### Phase 4: Documentation
- [ ] Update `config/config.yaml` comments
- [ ] Add benchmarking section to main README
- [ ] Create example truth sets for common organisms

---

## 9. API Reference

### 9.1 Helper Functions (common.smk)

```python
def get_all_result_tsvs(wildcards):
    """
    Collect all *_results.tsv paths for activated tools.
    Returns list of file paths for both per-comparison and per-sample tools.
    """
```

### 9.2 Script Functions (accuracy_benchmark.py)

```python
def tool_from_path(path: str) -> str:
    """Extract tool name from result file path."""

def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Map various column names to standard 'transcript' and 'position'."""

def detect_score_column(df: pd.DataFrame, tool_name: str = None) -> str:
    """Detect the most likely score column for ranking predictions."""

def compute_mcc(tp: int, fp: int, fn: int, tn: int) -> float:
    """Compute Matthews Correlation Coefficient."""

def compute_ranking_metrics(pred_df, truth_pos, truth_neg, score_col, window) -> tuple:
    """Compute AUPRC and AUROC given predictions with scores."""
```

---

## 10. Error Handling

### 10.1 Missing Truth Set
- If `benchmark.truth_set` is empty string → Skip accuracy benchmarking
- If file does not exist → Raise error in `accuracy_benchmark` rule

### 10.2 Missing Columns
- Truth set missing required columns → Python assertion error with clear message
- Tool output missing transcript/position → Warning logged, file skipped

### 10.3 Empty Results
- No predictions from any tool → Empty output file with header only
- No truth sites for modification type → Skipped for that type

### 10.4 Inferred Negatives
- When `label` column is missing, negatives are inferred from predictions
- `specificity` and `mcc` are set to `NaN` (not meaningful for inferred negatives)

---

## 11. Testing Strategy

### 11.1 Unit Tests
- Test `normalize_columns()` with various column name variants
- Test `detect_score_column()` with tool-specific dataframes
- Test `compute_mcc()` with edge cases (zeros, perfect scores)

### 11.2 Integration Tests
- Run full benchmarking with sample truth set in `.test/`
- Verify output file format and metric calculations
- Test with multiple window values

### 11.3 Test Data
- Synthetic truth set with known TP/FP/FN/TN counts
- Tool outputs with known score distributions
