# Per-Tool Benchmark Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor the monolithic `accuracy_benchmark` rule into per-tool parallel benchmark rules with coverage-first fair comparison.

**Architecture:** 5-layer DAG (coverage → union → native/fair eval → aggregation → downstream path updates) using `{tool}` and `{comparison}` Snakemake wildcards. Three new Python scripts replace the monolithic `accuracy_benchmark.py`. Downstream rules get path-only updates.

**Tech Stack:** Snakemake 6.4+, Python 3.10+, pandas, numpy, scikit-learn

**Spec:** `docs/superpowers/specs/2026-03-26-per-tool-benchmark-design.md`

---

### Task 1: Add PER_COMPARISON_TOOLS constant to common.smk

**Files:**
- Modify: `workflow/rules/common.smk:57` (after `comparisons` definition)
- Modify: `workflow/Snakefile:150-156` (wildcard_constraints)

- [ ] **Step 1: Add PER_COMPARISON_TOOLS list and helper**

In `workflow/rules/common.smk`, after the `comparisons` line (line 57), add:

```python
PER_COMPARISON_TOOLS = ["xpore", "nanocompore", "baleen", "differr", "drummer",
                        "eligos2", "epinano", "psipore", "pybaleen"]


def get_active_comparison_tools():
    return [t for t in config["tools"]
            if config["tools"][t]["activate"] and t in PER_COMPARISON_TOOLS]
```

- [ ] **Step 2: Add wildcard constraint for tool**

In `workflow/Snakefile`, add `tool` to the `wildcard_constraints` block (after line 156):

```python
    tool="[a-z0-9]+",
```

- [ ] **Step 3: Verify dry run still works**

Run: `snakemake --use-conda --dry-run --cores 1 2>&1 | tail -5`
Expected: No errors (existing rules unaffected by new constant)

- [ ] **Step 4: Commit**

```bash
git add workflow/rules/common.smk workflow/Snakefile
git commit -m "Add PER_COMPARISON_TOOLS constant and tool wildcard constraint"
```

---

### Task 2: Create benchmark_coverage.py script

**Files:**
- Create: `workflow/scripts/benchmark_coverage.py`

This script handles both per-tool coverage (Layer 1) and coverage union (Layer 2), detected by output keys.

- [ ] **Step 1: Create the script with per-tool coverage mode**

Create `workflow/scripts/benchmark_coverage.py` with two modes:

1. **Per-tool mode** (output has `covered` key): Load truth set and tool results, find which truth sites the tool can score within the positional window, output as `transcript\tposition` TSV.

2. **Union mode** (output has `union` and `called_sites` keys): Read all `{tool}_covered.tsv` files, take union of all (transcript, position) pairs, track which tools cover each site (alphabetical comma-separated `covered_by_tools` column), count covered sites per tool.

Key functions to inline from `benchmark_utils.py`:
- `normalize_columns()` — map variant column names to standard `transcript` and `position`
- `parse_windows()` — handle single int or list of windows from config

The per-tool coverage function builds a lookup of predicted positions by transcript, then for each truth site checks if any prediction exists within `max(windows)` nucleotides.

The union function reads all covered TSVs, collects (transcript, position) → set of tool names, writes union.tsv with `covered_by_tools` column and called_sites.tsv with per-tool counts.

Mode detection: check `snakemake.output` keys for `covered` vs `union`.

See spec Section 5.1 for complete logic.

- [ ] **Step 2: Verify syntax**

Run: `python -c "import py_compile; py_compile.compile('workflow/scripts/benchmark_coverage.py', doraise=True)"`
Expected: No output (clean compile)

- [ ] **Step 3: Commit**

```bash
git add workflow/scripts/benchmark_coverage.py
git commit -m "Add benchmark_coverage.py for per-tool and union coverage"
```

---

### Task 3: Create benchmark_per_tool.py script

**Files:**
- Create: `workflow/scripts/benchmark_per_tool.py`

Handles both native (Layer 3a) and fair (Layer 3b) evaluation. Detects mode by presence of `union` in `snakemake.input`.

- [ ] **Step 1: Create the script**

Create `workflow/scripts/benchmark_per_tool.py` with:

**Mode detection:** If `snakemake.input` has `union` key → fair mode, otherwise → native mode.

**Core logic (shared between modes):**
1. Load tool results and truth set, normalize columns
2. Detect candidate score columns using `TOOL_SCORE_MAP` (inline the full map from `accuracy_benchmark.py` lines 157-203, covering all 9 comparison-based tools)
3. For each candidate column that exists in the data:
   a. Extract scores for all evaluation sites
   b. In fair mode: sites in union but missing from tool predictions get score=0
   c. Detect if p-value column using `is_pvalue_column()` (inline from `benchmark_utils.py`)
   d. If p-value: apply -log10 transform (clamp to 1e-300 minimum)
   e. Compute AUROC and AUPRC via `sklearn.metrics.roc_auc_score` and `average_precision_score`
   f. Find optimal threshold maximizing F1 over 100 linearly-spaced thresholds
   g. Convert threshold back to original scale if p-value transformed
4. Rank all candidates by AUROC, pick best

**Outputs:**
- `score_comparison.tsv` — all candidates with columns: score_column, is_pvalue, transform, auroc, prauc, best_threshold, best_threshold_original, f1, precision, recall
- `best_metrics.tsv` — single row for the best column
- `best_score.txt` — plain text column name

If tool produces no valid results, write empty files with headers.

Key inline functions: `normalize_columns()`, `is_pvalue_column()`, `match_column()`, `tool_from_path()`, `parse_windows()`.

See spec Sections 4 and 5.2 for complete logic.

- [ ] **Step 2: Verify syntax**

Run: `python -c "import py_compile; py_compile.compile('workflow/scripts/benchmark_per_tool.py', doraise=True)"`
Expected: No output (clean compile)

- [ ] **Step 3: Commit**

```bash
git add workflow/scripts/benchmark_per_tool.py
git commit -m "Add benchmark_per_tool.py for native and fair evaluation"
```

---

### Task 4: Create benchmark_aggregation.py script

**Files:**
- Create: `workflow/scripts/benchmark_aggregation.py`

- [ ] **Step 1: Create the script**

Create `workflow/scripts/benchmark_aggregation.py` that:

1. Reads all `best_metrics.tsv` files from native and fair modes (extracting tool name and comparison from file path structure `.../native/{tool}/{comparison}/best_metrics.tsv`)
2. Uses native results as primary evaluation source
3. Produces these aggregated outputs matching downstream format expectations:

**accuracy_summary.tsv** — per tool × modification type. Since per-tool evaluation doesn't split by mod_type, use `modification_type='all'` with the tool's metrics from each comparison. Columns: tool, modification_type, window, precision, recall, f1, tp, fp, fn, tn, specificity, mcc, auprc, auroc, called_sites, total_truth, total_predicted, total_negative.

**accuracy_summary_overall.tsv** — per tool averaged across comparisons. Same columns minus modification_type.

**accuracy_summary_by_comparison.tsv** — per tool × comparison. Columns include comparison.

**accuracy_summary_by_negative_type.tsv** — empty DataFrame with correct columns (negative type analysis stays in its own rule).

**by_tool.tsv** (new) — mean±std metrics per tool: tool, best_score, mean_auroc, std_auroc, mean_prauc, std_prauc, mean_f1, std_f1, n_comparisons.

**best_scores.tsv** (new) — most common best score column per tool: tool, best_score, n_comparisons_best.

**called_sites_by_comparison.tsv** — aggregated from coverage called_sites.tsv files, adding comparison column.

**called_sites_summary.tsv** — per-tool totals from called_sites_by_comparison.

Touch `.benchmark_complete` marker for downstream rules.

See spec Section 5.3 for complete logic.

- [ ] **Step 2: Verify syntax**

Run: `python -c "import py_compile; py_compile.compile('workflow/scripts/benchmark_aggregation.py', doraise=True)"`
Expected: No output (clean compile)

- [ ] **Step 3: Commit**

```bash
git add workflow/scripts/benchmark_aggregation.py
git commit -m "Add benchmark_aggregation.py for cross-tool result aggregation"
```

---

### Task 5: Replace benchmark rules in benchmark_accuracy.smk

**Files:**
- Modify: `workflow/rules/benchmark_accuracy.smk` (full rewrite)

- [ ] **Step 1: Replace the file contents**

Replace `workflow/rules/benchmark_accuracy.smk` with the new per-tool rules. Keep `get_all_result_tsvs` (still used by other rules like kmer_negatives, multithreshold). Remove `accuracy_benchmark` rule and `get_benchmark_output_files` helper.

New rules (all guarded by `if config.get("benchmark", {}).get("truth_set", ""):`):

**benchmark_tool_coverage:**
- input: `{project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv`, truth_set
- output: `{project}/results/benchmarks/coverage/{comparison}/{tool}_covered.tsv`
- params: window from config
- priority: 30
- script: `../scripts/benchmark_coverage.py`

**benchmark_coverage_union:**
- input: lambda expanding all tool covered files for the comparison using `get_active_comparison_tools()`
- output: `{project}/results/benchmarks/coverage/{comparison}/union.tsv`, `called_sites.tsv`
- priority: 31
- script: `../scripts/benchmark_coverage.py`

**benchmark_tool_native:**
- input: tool results, truth_set (NO union dependency)
- output: `{project}/results/benchmarks/native/{tool}/{comparison}/score_comparison.tsv`, `best_metrics.tsv`, `best_score.txt`
- params: window
- priority: 32
- script: `../scripts/benchmark_per_tool.py`

**benchmark_tool_fair:**
- input: tool results, union.tsv, truth_set
- output: `{project}/results/benchmarks/fair/{tool}/{comparison}/score_comparison.tsv`, `best_metrics.tsv`, `best_score.txt`
- params: window
- priority: 33
- script: `../scripts/benchmark_per_tool.py`

**benchmark_aggregate:**
- input: lambda expanding all native and fair best_metrics.tsv, all called_sites.tsv, truth_set
- output: all aggregated files under `{project}/results/benchmarks/aggregated/`, plus `.benchmark_complete` touch file
- priority: 34
- script: `../scripts/benchmark_aggregation.py`

All rules include `log:`, `benchmark:`, `container:`, `resources:`, and `threads:` directives following existing patterns.

- [ ] **Step 2: Verify Snakemake can parse the file**

Run: `snakemake --lint --snakefile workflow/Snakefile 2>&1 | head -20`
Expected: No fatal errors

- [ ] **Step 3: Commit**

```bash
git add workflow/rules/benchmark_accuracy.smk
git commit -m "Replace monolithic accuracy_benchmark with per-tool benchmark rules"
```

---

### Task 6: Update downstream rule input paths

**Files:**
- Modify: `workflow/rules/benchmark_viz.smk`

- [ ] **Step 1: Update benchmark_statistics inputs**

Change inputs from `{project}/results/benchmarks/accuracy_summary.tsv` to `{project}/results/benchmarks/aggregated/accuracy_summary.tsv`. Same for `accuracy_summary_by_comparison.tsv`.

- [ ] **Step 2: Update benchmark_visualization inputs**

Remove the `benchmarks` directory input. Keep the `.benchmark_complete` touch file and resource_summary. Add explicit `summary="{project}/results/benchmarks/aggregated/accuracy_summary.tsv"` input.

- [ ] **Step 3: Update benchmark_r_figures inputs**

Change the `aggregated` expand from `benchmarks/{agg_file}` to `benchmarks/aggregated/{agg_file}`.

Change `optimal_scores` from `benchmarks/optimal_score_per_tool.tsv` to `benchmarks/aggregated/best_scores.tsv`.

Update `snakemake@input` key names in `run_all_figures.R` if they changed, or keep them compatible.

- [ ] **Step 4: Update benchmark_pdf_report inputs**

Remove the `benchmarks` directory input. Use explicit file inputs from `aggregated/` path. Update `params.benchmark_dir` to `{project}/results/benchmarks/aggregated`.

- [ ] **Step 5: Commit**

```bash
git add workflow/rules/benchmark_viz.smk
git commit -m "Update downstream benchmark rules to read from aggregated/ paths"
```

---

### Task 7: Update get_final_output() in common.smk

**Files:**
- Modify: `workflow/rules/common.smk:322-376` (benchmark output section)

- [ ] **Step 1: Replace the benchmark outputs block**

Replace the benchmark section in `get_final_output()` (starting at the comment `# Aggregated accuracy benchmarking outputs`) with new per-tool outputs:

1. Coverage files via `expand()` over comparisons × active_tools
2. Union files via `expand()` over comparisons
3. Native and fair best_metrics via `expand()` over tools × comparisons
4. Aggregated summaries under `aggregated/` path
5. Keep all unchanged downstream outputs (viz, reports, thresholds, negatives, statistics, sensitivity, figures) — only change paths that moved to `aggregated/`

Use `get_active_comparison_tools()` for tool lists.

- [ ] **Step 2: Verify dry run**

Run: `snakemake --use-conda --dry-run --cores 1 2>&1 | tail -10`
Expected: No errors about missing rules or circular dependencies

- [ ] **Step 3: Commit**

```bash
git add workflow/rules/common.smk
git commit -m "Update get_final_output() for per-tool benchmark structure"
```

---

### Task 8: Delete old accuracy_benchmark.py

**Files:**
- Delete: `workflow/scripts/accuracy_benchmark.py`

- [ ] **Step 1: Verify no remaining references**

Run: `grep -r "accuracy_benchmark" workflow/rules/ --include="*.smk"`
Expected: No references to `accuracy_benchmark` as a rule name or script reference

- [ ] **Step 2: Delete the old script**

```bash
git rm workflow/scripts/accuracy_benchmark.py
```

- [ ] **Step 3: Commit**

```bash
git commit -m "Remove old monolithic accuracy_benchmark.py"
```

---

### Task 9: Verify full pipeline with dry run

**Files:** None (verification only)

- [ ] **Step 1: Run dry run from project root**

Run: `snakemake --use-conda --dry-run --cores 1 2>&1 | grep -E "(benchmark_tool|benchmark_coverage|benchmark_aggregate)" | head -20`
Expected: See per-tool rules listed (benchmark_tool_coverage, benchmark_tool_native, benchmark_tool_fair, benchmark_coverage_union, benchmark_aggregate)

- [ ] **Step 2: Check rule count**

Run: `snakemake --use-conda --dry-run --cores 1 2>&1 | grep "rules" | tail -3`
Expected: Shows the total number of rules including the new benchmark rules

- [ ] **Step 3: Verify DAG structure**

Run: `snakemake --use-conda --dry-run --cores 1 --dag 2>/dev/null | head -50`
Expected: DAG output shows benchmark rules with correct dependencies (coverage → union → fair, native runs independently)

- [ ] **Step 4: Run test suite dry run**

Run: `cd .test && snakemake --use-conda --dry-run --cores 1 2>&1 | tail -10`
Expected: No errors. The test config has `benchmark.truth_set` set, so per-tool benchmark rules should appear in the DAG.

- [ ] **Step 5: Commit if fixes needed**

If any fixes were needed during verification:
```bash
git add -u
git commit -m "Fix issues found during dry-run verification"
```
