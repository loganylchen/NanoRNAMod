# NanoRNAMod R Visualization Scripts

Nature-quality ggplot2 visualization scripts for benchmark figures.

## Overview

This directory contains R scripts for generating publication-ready benchmark figures following Nature journal style guidelines. All figures use colorblind-friendly palettes (Okabe-Ito) and are output at 300 DPI resolution.

## Scripts

### Utility Script
- **`00_utils.R`** - Shared utilities, themes, and helper functions

### Figure Scripts
1. **`01_performance_overview.R`** - Main performance comparison
   - F1 score bar chart (tools sorted by performance)
   - Precision and Recall grouped bar chart
   - AUROC bar chart (if available)
   - AUPRC bar chart (if available)

2. **`02_per_modification_analysis.R`** - Stratified by modification type
   - Heatmap of F1 scores (tools × modification types)
   - Grouped bar chart per modification type
   - Faceted precision-recall plots
   - Window-based faceted plots

3. **`03_threshold_optimization.R`** - Threshold analysis
   - Precision-recall tradeoff curves
   - Faceted metric curves by tool
   - F1 vs threshold with optimal highlighted
   - Threshold stability analysis

4. **`04_resource_comparison.R`** - Computational resources
   - Runtime comparison (linear and log scale)
   - Memory usage comparison
   - I/O comparison
   - Efficiency scatter plots (F1 vs resources)

5. **`05_optimal_metrics_summary.R`** - Summary visualizations
   - Optimal configuration table
   - Score column usage summary
   - Tool performance ranking
   - Best tool per metric highlight

### Master Script
- **`run_all_figures.R`** - Generate all figures with one command

## Usage

### Generate All Figures

```bash
# From benchmark results directory
Rscript workflow/scripts/R/run_all_figures.R \
  --data-dir results/benchmark \
  --output results/figures

# With specific window filter
Rscript workflow/scripts/R/run_all_figures.R \
  --data-dir results/benchmark \
  --output results/figures \
  --window 10
```

### Generate Individual Figures

```bash
# Performance overview
Rscript workflow/scripts/R/01_performance_overview.R \
  --input results/benchmark/accuracy_summary.tsv \
  --output results/figures/01_performance_overview

# Per-modification analysis
Rscript workflow/scripts/R/02_per_modification_analysis.R \
  --input results/benchmark/accuracy_summary.tsv \
  --output results/figures/02_per_modification

# Threshold optimization
Rscript workflow/scripts/R/03_threshold_optimization.R \
  --input results/benchmark/threshold_evaluation.tsv \
  --output results/figures/03_threshold_optimization

# Resource comparison
Rscript workflow/scripts/R/04_resource_comparison.R \
  --resources results/benchmark/resource_summary.tsv \
  --accuracy results/benchmark/accuracy_summary_overall.tsv \
  --output results/figures/04_resource_comparison

# Optimal metrics summary
Rscript workflow/scripts/R/05_optimal_metrics_summary.R \
  --optimal results/benchmark/optimal_thresholds.tsv \
  --accuracy results/benchmark/accuracy_summary_overall.tsv \
  --output results/figures/05_optimal_summary
```

## Input Data Format

### accuracy_summary.tsv / accuracy_summary_overall.tsv
```
tool	modification_type	window	precision	recall	f1	tp	fp	fn	auprc	auroc	mcc
xpore	m6A	0	0.85	0.78	0.81	125	22	35	0.82	0.88	0.68
nanocompore	m6A	0	0.72	0.85	0.78	136	53	25	0.79	0.85	0.61
...
```

### resource_summary.tsv
```
tool	sample	s	max_rss	cpu_time	io_in	io_out
xpore	sample1	245.3	1024	485.2	1024	512
nanocompore	sample1	512.7	2048	1024.5	2048	1024
...
```

### threshold_evaluation.tsv
```
tool	score_column	threshold	precision	recall	f1
xpore	score	0.5	0.70	0.85	0.76
xpore	score	0.6	0.78	0.80	0.79
...
```

### optimal_thresholds.tsv
```
tool	score_column	threshold	f1	precision	recall
xpore	score	0.65	0.82	0.85	0.79
nanocompore	probability	0.45	0.78	0.75	0.81
...
```

## Output

Each script generates both PDF (vector format for publication) and PNG (raster format for preview) versions at 300 DPI.

```
results/figures/
├── 01_performance_overview/
│   ├── 01_f1_scores.pdf
│   ├── 01_f1_scores.png
│   ├── 02_precision_recall.pdf
│   └── ...
├── 02_per_modification/
│   ├── modification_01_heatmap.pdf
│   └── ...
├── 03_threshold_optimization/
│   ├── threshold_01_pr_tradeoff.pdf
│   └── ...
├── 04_resource_comparison/
│   ├── resource_01_runtime.pdf
│   └── ...
└── 05_optimal_summary/
    ├── summary_01_optimal_table.pdf
    └── ...
```

## Style Guidelines

- **Font**: Helvetica (or Arial as fallback)
- **Base size**: 8pt for Nature single-column figures
- **Colors**: Okabe-Ito colorblind-friendly palette
- **Dimensions**:
  - Single column: 85mm width
  - Double column: 174mm width
  - Heights vary by content (70-120mm typical)
- **Resolution**: 300 DPI for PNG output

## Dependencies

Required R packages:
- ggplot2
- dplyr
- tidyr
- readr
- scales
- grid
- gridExtra
- RColorBrewer
- viridis

Install missing packages automatically by running the utility script:
```r
source("workflow/scripts/R/00_utils.R")
```

## Customization

To modify figure appearance, edit the `theme_nature()` function in `00_utils.R`:

```r
theme_nature <- function(base_size = 8, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.25),
      # ... other theme settings
    )
}
```

## Troubleshooting

### Font Not Available
If Helvetica is not available, the script will automatically fall back to Arial or sans-serif.

### Missing Data Columns
Scripts check for column availability and skip figures when data is missing. Warnings are printed for missing optional data.

### Memory Issues with Large Datasets
For very large threshold sweep files, consider:
1. Filtering to specific tools: `df <- df[df$tool %in% c("xpore", "nanocompore"), ]`
2. Aggregating before plotting
3. Using data.table instead of dplyr for large datasets
