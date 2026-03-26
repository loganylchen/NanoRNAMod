# Per-Tool Benchmark System Design (Revised)

**Date**: 2026-03-26 (revised)
**Status**: Approved
**Author**: Claude (with user collaboration)

## 1. Overview

### 1.1 Problem Statement

The current benchmark system uses a single aggregated rule (`accuracy_benchmark`) that processes all modification detection tools together. This design has limitations:

- **No parallelization**: All tools must be processed sequentially in one rule
- **No incremental updates**: Adding a new tool requires re-running the entire benchmark
- **Poor modularity**: Tool-specific logic is mixed with aggregation logic

### 1.2 Proposed Solution

Refactor into per-tool benchmark rules using `{tool}` wildcards that run independently and in parallel, with a coverage-first approach for fair cross-tool comparison.

### 1.3 Design Goals

- **Parallelization**: Each tool's benchmark runs independently via Snakemake wildcards
- **Incremental updates**: Re-run only affected tools when data changes
- **Modularity**: Clear separation between coverage, evaluation, and aggregation
- **Fair comparison**: Coverage union ensures all tools evaluated on same sites
- **Native evaluation**: Tools also evaluated on only their covered sites
- **Clean break**: Replaces old monolithic rule directly; no legacy adapters

### 1.4 Scope

**In scope**: 9 comparison-based tools (xpore, nanocompore, baleen, differr, drummer, eligos2, epinano, psipore, pybaleen)

**Out of scope**: Per-sample tools (tandemmod, directrm, m6atm, rnano, nanopsu, nanomud, penguin) — deferred to a future phase

## 2. Coverage-First Design

### 2.1 Core Concept

For each comparison (native_sample vs control_sample):

1. **Compute coverage**: Each tool identifies which truth sites it can score
2. **Union coverage**: Take the union of all tools' covered sites
3. **Native evaluation**: Each tool evaluated only on sites it covers
4. **Fair evaluation**: All tools evaluated on union sites; missing predictions treated as negative (score = 0)

### 2.2 Coverage Definition

A site is "covered" by a tool if the tool produces any prediction (with score) for that transcript/position within the configured positional tolerance window. No distinction between modification types.

### 2.3 Evaluation Modes

| Mode | Evaluation Sites | Missing Predictions |
|------|-----------------|---------------------|
| Native | Tool's covered sites only | Excluded from calculation |
| Fair | Union of all tools' covered sites | Treated as negative (score = 0) |

## 3. Rule Structure

### 3.1 Execution DAG

```
Layer 1: benchmark_tool_coverage        (per tool × per comparison, parallel)
Layer 2: benchmark_coverage_union        (per comparison, waits for all L1)
Layer 3a: benchmark_tool_native          (per tool × per comparison, parallel)
Layer 3b: benchmark_tool_fair            (per tool × per comparison, parallel, depends on L2)
Layer 4: benchmark_aggregate             (single rule, collects all L3 outputs)
Layer 5: existing downstream             (statistics, sensitivity, R figures — path updates only)
```

All rules use `{tool}` and `{comparison}` wildcards. Snakemake resolves via `expand()` in `get_final_output()`.

`PER_COMPARISON_TOOLS` is defined in `common.smk`:
```python
PER_COMPARISON_TOOLS = ["xpore", "nanocompore", "baleen", "differr", "drummer",
                        "eligos2", "epinano", "psipore", "pybaleen"]
```

### 3.2 Layer 1: benchmark_tool_coverage

```python
rule benchmark_tool_coverage:
    input:
        results="{project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv",
        truth_set=config["benchmark"]["truth_set"],
    output:
        covered="{project}/results/benchmarks/coverage/{comparison}/{tool}_covered.tsv",
    params:
        window=config["benchmark"]["window"],
    resources:
        mem_mb=1024 * 4,
    threads: 1
    priority: 30
    log:
        "logs/{project}/benchmark_tool_coverage/{comparison}/{tool}.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_coverage.py"
```

**Output format** (`{tool}_covered.tsv`):
```
transcript	position
ENST0001	1234
ENST0001	5678
```

**Logic**:
1. Load truth set
2. Load tool predictions
3. For each truth site, check if tool has prediction within window
4. Output covered sites

### 3.3 Layer 2: benchmark_coverage_union

```python
rule benchmark_coverage_union:
    input:
        covered=lambda wc: expand(
            "{project}/results/benchmarks/coverage/{comparison}/{tool}_covered.tsv",
            project=wc.project, comparison=wc.comparison,
            tool=[t for t in config["tools"] if config["tools"][t]["activate"]
                  and t in PER_COMPARISON_TOOLS],
        ),
    output:
        union="{project}/results/benchmarks/coverage/{comparison}/union.tsv",
        called_sites="{project}/results/benchmarks/coverage/{comparison}/called_sites.tsv",
    resources:
        mem_mb=1024 * 2,
    threads: 1
    priority: 31
    log:
        "logs/{project}/benchmark_coverage_union/{comparison}.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_coverage.py"
```

**Output format** (`union.tsv`):
```
transcript	position	covered_by_tools
ENST0001	1234	baleen,nanocompore,xpore
ENST0001	5678	baleen,xpore
```

`covered_by_tools` is comma-separated, alphabetical order for reproducibility.

**Output format** (`called_sites.tsv`):
```
tool	n_covered
baleen	142
nanocompore	98
xpore	115
```

**Logic**:
1. Collect all covered sites from each tool's `_covered.tsv`
2. Union across all tools
3. Track which tools cover each site
4. Count covered sites per tool

**Mode detection**: The script detects union mode by checking for `union` in `snakemake.output` keys.

### 3.4 Layer 3a: benchmark_tool_native

```python
rule benchmark_tool_native:
    input:
        results="{project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv",
        truth_set=config["benchmark"]["truth_set"],
    output:
        scores="{project}/results/benchmarks/native/{tool}/{comparison}/score_comparison.tsv",
        best_metrics="{project}/results/benchmarks/native/{tool}/{comparison}/best_metrics.tsv",
        best_score="{project}/results/benchmarks/native/{tool}/{comparison}/best_score.txt",
    params:
        window=config["benchmark"]["window"],
    resources:
        mem_mb=1024 * 8,
    threads: 1
    priority: 32
    log:
        "logs/{project}/benchmark_tool_native/{tool}/{comparison}.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_per_tool.py"
```

**Can run in parallel with Layer 2** — does not depend on coverage union.

### 3.5 Layer 3b: benchmark_tool_fair

```python
rule benchmark_tool_fair:
    input:
        results="{project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv",
        union="{project}/results/benchmarks/coverage/{comparison}/union.tsv",
        truth_set=config["benchmark"]["truth_set"],
    output:
        scores="{project}/results/benchmarks/fair/{tool}/{comparison}/score_comparison.tsv",
        best_metrics="{project}/results/benchmarks/fair/{tool}/{comparison}/best_metrics.tsv",
        best_score="{project}/results/benchmarks/fair/{tool}/{comparison}/best_score.txt",
    params:
        window=config["benchmark"]["window"],
    resources:
        mem_mb=1024 * 8,
    threads: 1
    priority: 33
    log:
        "logs/{project}/benchmark_tool_fair/{tool}/{comparison}.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_per_tool.py"
```

**Depends on Layer 2**: Must wait for coverage union.

**Mode detection**: Script checks if `union` is present in `snakemake.input` keys. If present → fair mode (load union, fill missing with score=0). If absent → native mode (evaluate tool's own sites only).

### 3.6 Layer 4: benchmark_aggregate

```python
rule benchmark_aggregate:
    input:
        native=lambda wc: expand(
            "{project}/results/benchmarks/native/{tool}/{comparison}/best_metrics.tsv",
            project=wc.project,
            tool=[t for t in config["tools"] if config["tools"][t]["activate"]
                  and t in PER_COMPARISON_TOOLS],
            comparison=comparisons,
        ),
        fair=lambda wc: expand(
            "{project}/results/benchmarks/fair/{tool}/{comparison}/best_metrics.tsv",
            project=wc.project,
            tool=[t for t in config["tools"] if config["tools"][t]["activate"]
                  and t in PER_COMPARISON_TOOLS],
            comparison=comparisons,
        ),
        called_sites=lambda wc: expand(
            "{project}/results/benchmarks/coverage/{comparison}/called_sites.tsv",
            project=wc.project,
            comparison=comparisons,
        ),
        truth_set=config["benchmark"]["truth_set"],
    output:
        summary="{project}/results/benchmarks/aggregated/accuracy_summary.tsv",
        overall="{project}/results/benchmarks/aggregated/accuracy_summary_overall.tsv",
        by_comparison="{project}/results/benchmarks/aggregated/accuracy_summary_by_comparison.tsv",
        by_negative_type="{project}/results/benchmarks/aggregated/accuracy_summary_by_negative_type.tsv",
        by_tool="{project}/results/benchmarks/aggregated/by_tool.tsv",
        best_scores="{project}/results/benchmarks/aggregated/best_scores.tsv",
        called_sites_comp="{project}/results/benchmarks/aggregated/called_sites_by_comparison.tsv",
        called_sites_sum="{project}/results/benchmarks/aggregated/called_sites_summary.tsv",
        done=touch("{project}/results/benchmarks/.benchmark_complete"),
    resources:
        mem_mb=1024 * 4,
    threads: 1
    priority: 34
    log:
        "logs/{project}/benchmark_aggregate/aggregate.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_aggregation.py"
```

## 4. Output Formats

### 4.1 score_comparison.tsv (per-tool, per-comparison)

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

### 4.2 best_metrics.tsv (per-tool, per-comparison)

Single row with the best score column's metrics. Same columns as score_comparison.tsv.

### 4.3 best_score.txt (per-tool, per-comparison)

Plain text file containing only the score column name that achieved highest AUROC. No header.

### 4.4 Aggregated Outputs

`accuracy_summary.tsv` — per tool × modification type, matching the format consumed by downstream rules (statistics, sensitivity, R figures).

`accuracy_summary_overall.tsv` — per tool, averaged across comparisons and modification types.

`accuracy_summary_by_comparison.tsv` — per tool × comparison.

`accuracy_summary_by_negative_type.tsv` — per negative control strategy.

`by_tool.tsv`:
```
tool	best_score	mean_auroc	std_auroc	mean_prauc	std_prauc	mean_f1	std_f1	n_comparisons
xpore	p_value	0.89	0.05	0.85	0.06	0.82	0.04	3
```

`best_scores.tsv`:
```
tool	best_score	n_comparisons_best
xpore	p_value	2
nanocompore	logit_pvalue	1
```

## 5. Script Design

### 5.1 benchmark_coverage.py

Two modes, detected by output keys:

**Per-tool mode** (output has `covered`):
1. Load truth set, load tool results TSV
2. Normalize columns using inlined helpers from benchmark_utils.py
3. For each truth site, check if tool has a prediction within window nucleotides
4. Output matched sites as `transcript\tposition` TSV

**Union mode** (output has `union` and `called_sites`):
1. Read all `{tool}_covered.tsv` files from the input list
2. Union all (transcript, position) pairs
3. Track which tools cover each site → `covered_by_tools` column (alphabetical)
4. Write `union.tsv`
5. Count covered sites per tool → `called_sites.tsv`

### 5.2 benchmark_per_tool.py

Single script, two modes detected by presence of `union` in inputs:

**Core logic (shared)**:
1. Load tool results, load truth set
2. Detect candidate score columns using `tool_score_map`
3. For each score column: type conversion, p-value detection, -log10 transform
4. Compute AUROC, AUPRC via sklearn
5. Find optimal threshold (maximize F1) with 100 linearly-spaced thresholds
6. Compute precision/recall/F1 at optimal threshold
7. Rank columns by AUROC, pick best
8. Write score_comparison.tsv, best_metrics.tsv, best_score.txt

**Native mode** (no `union` in input):
- Evaluate only on sites where tool has predictions matching truth sites within window

**Fair mode** (`union` in input):
- Load union sites
- For union sites not in tool's predictions: assign score = 0
- Evaluate on all union sites

**P-value handling**:
```python
def is_pvalue_column(col_name):
    pvalue_keywords = ['p_value', 'pvalue', 'pval', 'p.value']
    col_lower = col_name.lower().replace(' ', '_')
    return any(kw in col_lower for kw in pvalue_keywords)
```

P-value columns get -log10 transformed so higher = better. Thresholds are converted back to original scale for reporting.

**Score column candidates** (from existing tool_score_map):
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
    'pybaleen': ['mod_score', 'score'],
}
```

Only columns that actually exist in the tool's output are evaluated. Missing columns are silently skipped.

### 5.3 benchmark_aggregation.py

1. Collect all `best_metrics.tsv` from native and fair modes across tools and comparisons
2. Parse tool name and comparison from file paths
3. Produce aggregated summaries matching the format downstream rules expect
4. Aggregate called_sites from coverage layer
5. Touch `.benchmark_complete` marker

**Import strategy**: All three scripts inline the needed utility functions (same pattern used by existing benchmark scripts to avoid Singularity/Snakemake import issues).

## 6. Directory Structure

```
{project}/results/benchmarks/
├── coverage/
│   └── {comparison}/
│       ├── {tool}_covered.tsv
│       ├── union.tsv
│       └── called_sites.tsv
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
│   ├── accuracy_summary.tsv
│   ├── accuracy_summary_overall.tsv
│   ├── accuracy_summary_by_comparison.tsv
│   ├── accuracy_summary_by_negative_type.tsv
│   ├── by_tool.tsv
│   ├── best_scores.tsv
│   ├── called_sites_by_comparison.tsv
│   └── called_sites_summary.tsv
├── .benchmark_complete
├── statistics/                        (unchanged)
├── sensitivity/                       (unchanged)
├── figures/                           (unchanged)
├── data/                              (unchanged)
└── viz/                               (unchanged)
```

## 7. Downstream Rule Updates

Path-only changes — no logic modifications.

| Rule | Change |
|------|--------|
| `benchmark_statistics` | Input paths: `benchmarks/` → `benchmarks/aggregated/` |
| `benchmark_sensitivity` | No change (reads `detailed_predictions.tsv` from `benchmark_detailed` rule) |
| `benchmark_visualization` | Input: explicit aggregated file paths instead of directory |
| `benchmark_r_figures` | Input paths: `benchmarks/` → `benchmarks/aggregated/` |
| `benchmark_pdf_report` | Input: explicit aggregated file paths instead of directory |
| `benchmark_kmer_negatives` | No change (reads raw tool results directly) |
| `benchmark_same_base_negatives` | No change |
| `benchmark_multithreshold` | No change |
| `benchmark_score_optimization` | No change |
| `benchmark_threshold` | No change |
| `benchmark_detailed` | No change |

### get_final_output() changes

Replace current benchmark output list with:
- Per-tool coverage files via `expand()`
- Per-comparison union and called_sites files
- Per-tool native/fair metrics via `expand()`
- Aggregated summary files under `aggregated/`
- All downstream outputs (statistics, sensitivity, figures) with updated paths where needed

### Deletions

- Old `accuracy_benchmark` rule in `benchmark_accuracy.smk`
- `get_benchmark_output_files()` helper function
- Old per-tool subdirectory outputs pattern (`benchmarks/{tool}/accuracy_summary.tsv`)

## 8. Wildcard Constraints

Add to Snakefile wildcard_constraints:
```python
wildcard_constraints:
    tool="[a-z0-9]+",
```

This prevents the `{tool}` wildcard from matching path separators or unexpected patterns in benchmark rules.

## 9. Configuration

No config schema changes needed. The existing `benchmark` section already has all required fields:
```yaml
benchmark:
  truth_set: "config/benchmark_ecoli_rRNA.tsv"
  window: [0, 1, 2, 5]
  n_thresholds: 50
  custom_thresholds: []
```

## 10. Testing Strategy

### Integration test
- Run full pipeline on `.test/` dataset
- Verify all output files created with correct format
- Verify parallelization works (multiple tools execute concurrently)

### Validation
- Compare aggregated output metrics with old `accuracy_benchmark` results on same data
- Verify AUROC/AUPRC calculations match sklearn
- Confirm threshold conversion (transformed ↔ original scale) is correct
- Check that fair mode score=0 assignment produces expected penalty behavior
