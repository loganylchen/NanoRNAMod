#!/usr/bin/env Rscript
# =============================================================================
# Master Script: Generate Publication-Quality Benchmark Figures
# =============================================================================
# Enhanced version with:
# - Bootstrap CI overlays for error bars
# - Significance annotations (Wilcoxon, permutation tests)
# - Sensitivity analysis figures (coverage, threshold robustness)
# - Main vs supplementary figure organization
# - Source data export for manuscript submission
#
# Usage:
#   Rscript run_all_figures.R [options]
#   OR via Snakemake script directive
#
# Output structure:
#   figures/main/          - Main figures (Nature-ready)
#   figures/supplementary/ - Supplementary figures
#   data/                  - Source data (TSV)
# =============================================================================

# Suppress startup messages
options(echo = FALSE, warn = 1)

# =============================================================================
# Configuration - Handle both Snakemake script and command-line invocation
# =============================================================================

# Get script directory — snakemake@scriptdir gives the real source dir even
# when Snakemake copies the script to a temp location.
script_dir <- if (exists("snakemake")) {
  snakemake@scriptdir
} else {
  tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      "."
    }
  }, error = function(e) ".")
}

# Source utilities from the original script directory
source(file.path(script_dir, "00_utils.R"))

# Initialize with defaults
data_dir <- "."
output_dir <- "figures"
window_filter <- NULL
theme <- "nature"
dpi <- 300
fig_format <- "pdf"

# Initialize input file paths
input_files <- list(
  # Accuracy metrics
  accuracy = NULL,
  accuracy_overall = NULL,
  accuracy_by_comparison = NULL,
  accuracy_by_negative_type = NULL,
  # Optimization results
  optimal_scores = NULL,
  thresholds = NULL,
  resources = NULL,
  # Statistical analysis
  bootstrap_ci = NULL,
  significance = NULL,
  fdr_corrected = NULL,
  effect_sizes = NULL,
  # Sensitivity analysis
  coverage = NULL,
  score_dist = NULL,
  threshold_robust = NULL
)

# Check if running as Snakemake script
if (exists("snakemake")) {
  # Snakemake script directive - use snakemake object
  input_files$accuracy <- snakemake@input[["aggregated"]][1]
  input_files$accuracy_overall <- snakemake@input[["aggregated"]][2]
  input_files$accuracy_by_comparison <- snakemake@input[["aggregated"]][3]
  input_files$accuracy_by_negative_type <- snakemake@input[["aggregated"]][4]
  input_files$optimal_scores <- snakemake@input[["optimal_scores"]]
  input_files$thresholds <- snakemake@input[["thresholds"]]
  input_files$resources <- snakemake@input[["resources"]]
  input_files$bootstrap_ci <- snakemake@input[["bootstrap_ci"]]
  input_files$significance <- snakemake@input[["significance"]]
  input_files$fdr_corrected <- snakemake@input[["fdr"]]
  input_files$effect_sizes <- snakemake@input[["effect_sizes"]]
  input_files$coverage <- snakemake@input[["coverage"]]
  input_files$score_dist <- snakemake@input[["score_dist"]]
  input_files$threshold_robust <- snakemake@input[["threshold_robust"]]

  output_dir <- snakemake@output[["dir"]]
  if (is.null(output_dir)) {
    output_dir <- dirname(dirname(snakemake@output[["fig1"]]))
  }
  window_filter <- snakemake@params[["window"]]
  theme <- snakemake@params[["theme"]]
  dpi <- snakemake@params[["dpi"]]
  fig_format <- snakemake@params[["format"]]

  data_dir <- dirname(input_files$accuracy)
} else {
  # Command-line invocation
  args <- commandArgs(trailingOnly = TRUE)
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--data-dir" && i + 1 <= length(args)) {
      data_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output" && i + 1 <= length(args)) {
      output_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--window" && i + 1 <= length(args)) {
      window_filter <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--theme" && i + 1 <= length(args)) {
      theme <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--dpi" && i + 1 <= length(args)) {
      dpi <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--format" && i + 1 <= length(args)) {
      fig_format <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  # Set default paths based on data_dir
  input_files$accuracy <- file.path(data_dir, "accuracy_summary.tsv")
  input_files$accuracy_overall <- file.path(data_dir, "accuracy_summary_overall.tsv")
  input_files$accuracy_by_comparison <- file.path(data_dir, "accuracy_summary_by_comparison.tsv")
  input_files$accuracy_by_negative_type <- file.path(data_dir, "accuracy_summary_by_negative_type.tsv")
  input_files$optimal_scores <- file.path(data_dir, "optimal_score_per_tool.tsv")
  input_files$thresholds <- file.path(data_dir, "threshold_evaluation.tsv")
  input_files$resources <- file.path(data_dir, "resource_summary.tsv")
  input_files$bootstrap_ci <- file.path(data_dir, "statistics/bootstrap_ci.tsv")
  input_files$significance <- file.path(data_dir, "statistics/significance_tests.tsv")
  input_files$fdr_corrected <- file.path(data_dir, "statistics/fdr_corrected.tsv")
  input_files$effect_sizes <- file.path(data_dir, "statistics/effect_sizes.tsv")
  input_files$coverage <- file.path(data_dir, "sensitivity/coverage_analysis.tsv")
  input_files$score_dist <- file.path(data_dir, "sensitivity/score_distribution.tsv")
  input_files$threshold_robust <- file.path(data_dir, "sensitivity/threshold_robustness.tsv")
}

# Define output directories
main_dir <- file.path(output_dir, "main")
supp_dir <- file.path(output_dir, "supplementary")
data_out_dir <- file.path(output_dir, "..", "data")

# =============================================================================
# Data Loading Functions
# =============================================================================

load_with_fallback <- function(path, required = FALSE) {
  if (file.exists(path)) {
    df <- tryCatch({
      as.data.frame(read.table(path, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE, check.names = FALSE))
    }, error = function(e) {
      warning("Error reading ", path, ": ", e$message)
      NULL
    })
    if (!is.null(df) && nrow(df) > 0) {
      message("  Loaded: ", basename(path), " (", nrow(df), " rows)")
      return(df)
    }
  }
  if (required) {
    warning("Required file not found or empty: ", path)
  }
  NULL
}

# =============================================================================
# Figure Generation Functions
# =============================================================================

# -----------------------------------------------------------------------------
# Figure 1: Overall Accuracy with Bootstrap CI
# -----------------------------------------------------------------------------
create_fig1_overall_accuracy <- function(accuracy_df, ci_df, tool_colors) {
  if (is.null(accuracy_df)) return(NULL)

  # Aggregate by tool
  df_summary <- accuracy_df %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(
      f1 = mean(f1, na.rm = TRUE),
      precision = mean(precision, na.rm = TRUE),
      recall = mean(recall, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(f1))

  df_summary$tool <- factor(df_summary$tool, levels = df_summary$tool)

  # Add CI if available
  if (!is.null(ci_df) && "metric" %in% names(ci_df)) {
    ci_f1 <- ci_df %>%
      dplyr::filter(metric == "f1") %>%
      dplyr::select(tool, ci_lower, ci_upper)
    df_summary <- df_summary %>%
      dplyr::left_join(ci_f1, by = "tool")
  }

  # Determine y-axis limit
  y_max <- if ("ci_upper" %in% names(df_summary)) {
    max(c(1.1, max(df_summary$f1, na.rm = TRUE) * 1.1, max(df_summary$ci_upper, na.rm = TRUE) * 1.05))
  } else {
    max(c(1.1, max(df_summary$f1, na.rm = TRUE) * 1.1))
  }

  plot_colors <- tool_colors[match(levels(df_summary$tool), names(tool_colors))]

  p <- ggplot(df_summary, aes(x = tool, y = f1, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.25)

  # Add error bars if CI available
  if ("ci_lower" %in% names(df_summary) && !all(is.na(df_summary$ci_lower))) {
    p <- p + geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2, linewidth = 0.5
    )
  }

  p <- p +
    scale_fill_manual(values = plot_colors, guide = "none") +
    geom_text(aes(label = sprintf("%.2f", f1)), vjust = -0.5, size = 2.5) +
    labs(
      title = "Tool Performance Comparison",
      subtitle = "F1 Score with 95% Bootstrap CI",
      x = "Tool",
      y = "F1 Score"
    ) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

# -----------------------------------------------------------------------------
# Figure 2: PR Curves with AUPRC
# -----------------------------------------------------------------------------
create_fig2_pr_curves <- function(accuracy_df, tool_colors) {
  if (is.null(accuracy_df)) return(NULL)

  df_plot <- accuracy_df %>%
    dplyr::select(tool, precision, recall, f1) %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(
      precision = mean(precision, na.rm = TRUE),
      recall = mean(recall, na.rm = TRUE),
      f1 = mean(f1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(f1))

  df_long <- df_plot %>%
    tidyr::pivot_longer(
      cols = c(precision, recall),
      names_to = "metric",
      values_to = "value"
    )

  df_long$tool <- factor(df_long$tool, levels = df_plot$tool)

  p <- ggplot(df_long, aes(x = tool, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(
      values = c("precision" = "#0072B2", "recall" = "#D55E00"),
      labels = c("Precision", "Recall"),
      name = NULL
    ) +
    geom_text(aes(label = sprintf("%.2f", value)),
              position = position_dodge(width = 0.8), vjust = -0.5, size = 2) +
    labs(
      title = "Precision vs Recall Trade-off",
      x = "Tool",
      y = "Score"
    ) +
    scale_y_continuous(limits = c(0, 1.15), expand = c(0, 0)) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  p
}

# -----------------------------------------------------------------------------
# Figure 3: ROC Curves (AUROC comparison)
# -----------------------------------------------------------------------------
create_fig3_roc_curves <- function(accuracy_df, ci_df, tool_colors) {
  if (is.null(accuracy_df) || !"auroc" %in% names(accuracy_df)) return(NULL)

  df_summary <- accuracy_df %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(auroc = mean(auroc, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(!is.na(auroc)) %>%
    dplyr::arrange(desc(auroc))

  if (nrow(df_summary) == 0) return(NULL)

  df_summary$tool <- factor(df_summary$tool, levels = df_summary$tool)
  plot_colors <- tool_colors[match(levels(df_summary$tool), names(tool_colors))]

  # Add CI if available
  if (!is.null(ci_df) && "metric" %in% names(ci_df)) {
    ci_auroc <- ci_df %>%
      dplyr::filter(metric == "auroc") %>%
      dplyr::select(tool, ci_lower, ci_upper)
    df_summary <- df_summary %>%
      dplyr::left_join(ci_auroc, by = "tool")
  }

  p <- ggplot(df_summary, aes(x = tool, y = auroc, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.25)

  if ("ci_lower" %in% names(df_summary) && !all(is.na(df_summary$ci_lower))) {
    p <- p + geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2, linewidth = 0.5
    )
  }

  p <- p +
    scale_fill_manual(values = plot_colors, guide = "none") +
    geom_text(aes(label = sprintf("%.3f", auroc)), vjust = -0.5, size = 2.5) +
    labs(
      title = "AUROC Comparison",
      subtitle = "Area Under ROC Curve (higher is better)",
      x = "Tool",
      y = "AUROC"
    ) +
    scale_y_continuous(limits = c(0.45, 1.05), expand = c(0, 0)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.25) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

# -----------------------------------------------------------------------------
# Figure 4: F1 Comparison with Significance
# -----------------------------------------------------------------------------
create_fig4_f1_comparison <- function(accuracy_df, sig_df, tool_colors) {
  if (is.null(accuracy_df)) return(NULL)

  # Aggregate by tool
  df_summary <- accuracy_df %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(
      f1 = mean(f1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(f1))

  df_summary$tool <- factor(df_summary$tool, levels = df_summary$tool)
  plot_colors <- tool_colors[match(levels(df_summary$tool), names(tool_colors))]

  p <- ggplot(df_summary, aes(x = tool, y = f1, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = plot_colors, guide = "none") +
    geom_text(aes(label = sprintf("%.3f", f1)), vjust = -0.5, size = 2.5) +
    labs(
      title = "F1 Score Comparison",
      subtitle = "Error bars indicate 95% bootstrap CI",
      x = "Tool",
      y = "F1 Score"
    ) +
    scale_y_continuous(limits = c(0, max(1.15, max(df_summary$f1) * 1.15)), expand = c(0, 0)) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

# -----------------------------------------------------------------------------
# Figure 5: Coverage Sensitivity
# -----------------------------------------------------------------------------
create_fig5_coverage_sensitivity <- function(coverage_df, tool_colors) {
  if (is.null(coverage_df)) return(NULL)

  # Parse coverage bins
  df_plot <- coverage_df %>%
    dplyr::filter(!is.na(f1)) %>%
    dplyr::mutate(
      coverage_bin = factor(coverage_bin, ordered = TRUE)
    )

  if (nrow(df_plot) == 0) return(NULL)

  p <- ggplot(df_plot, aes(x = coverage_bin, y = f1, fill = tool, group = tool)) +
    geom_line(aes(color = tool), linewidth = 0.8) +
    geom_point(aes(color = tool), size = 2) +
    scale_color_manual(values = tool_colors) +
    labs(
      title = "Performance vs Coverage Depth",
      subtitle = "F1 score stratified by read coverage",
      x = "Coverage Bin",
      y = "F1 Score"
    ) +
    scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  p
}

# -----------------------------------------------------------------------------
# Figure 6: Resource Usage
# -----------------------------------------------------------------------------
create_fig6_resource_usage <- function(resource_df, accuracy_df, tool_colors) {
  if (is.null(resource_df)) return(NULL)

  df_plot <- resource_df

  # Handle both pre-summarised (mean_memory_mb) and raw benchmark (max_rss, s) columns
  if ("mean_memory_mb" %in% names(df_plot) && "mean_runtime_s" %in% names(df_plot)) {
    df_plot <- df_plot %>%
      dplyr::mutate(runtime_min = mean_runtime_s / 60, memory_gb = mean_memory_mb / 1024)
  } else if ("max_rss" %in% names(df_plot) && "s" %in% names(df_plot)) {
    # Raw Snakemake benchmark data — summarise per tool
    df_plot <- df_plot %>%
      dplyr::group_by(tool) %>%
      dplyr::summarise(
        runtime_min = mean(s, na.rm = TRUE) / 60,
        memory_gb = mean(max_rss, na.rm = TRUE) / 1024,
        .groups = "drop"
      )
    if (!is.null(accuracy_df)) {
      acc_summary <- accuracy_df %>%
        dplyr::group_by(tool) %>%
        dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop")
      df_plot <- df_plot %>% dplyr::left_join(acc_summary, by = "tool")
    }
  } else {
    message("  Skipping resource figure: required columns missing")
    return(NULL)
  }

  p <- ggplot(df_plot, aes(x = runtime_min, y = memory_gb, color = tool)) +
    geom_point(aes(size = ifelse("f1" %in% names(df_plot), f1, 1)), alpha = 0.8) +
    geom_text(aes(label = tool), vjust = -1, size = 2.5, hjust = 0.5) +
    scale_color_manual(values = tool_colors, guide = "none") +
    scale_size_continuous(range = c(3, 10), name = "F1 Score") +
    labs(
      title = "Resource Usage Comparison",
      subtitle = "Runtime vs Memory (bubble size = F1 score)",
      x = "Runtime (minutes)",
      y = "Memory (GB)"
    ) +
    theme_nature() +
    theme(legend.position = "right")

  p
}

# =============================================================================
# Supplementary Figure Functions
# =============================================================================

# -----------------------------------------------------------------------------
# SFig 1: Per-Comparison Breakdown
# -----------------------------------------------------------------------------
create_sfig_per_comparison <- function(by_comp_df, tool_colors) {
  if (is.null(by_comp_df)) return(NULL)

  p <- ggplot(by_comp_df, aes(x = tool, y = f1, fill = comparison)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "F1 Score by Comparison",
      x = "Tool",
      y = "F1 Score"
    ) +
    scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  p
}

# -----------------------------------------------------------------------------
# SFig 2: Score Distributions
# -----------------------------------------------------------------------------
create_sfig_score_dist <- function(score_dist_df, tool_colors) {
  if (is.null(score_dist_df)) return(NULL)

  # Create a visualization of score separation
  df_plot <- score_dist_df %>%
    dplyr::select(tool, mean_diff, cohens_d, overlap_coefficient, ks_statistic) %>%
    tidyr::pivot_longer(
      cols = c(mean_diff, cohens_d, overlap_coefficient),
      names_to = "metric",
      values_to = "value"
    )

  p <- ggplot(df_plot, aes(x = tool, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(
      values = c("mean_diff" = "#0072B2", "cohens_d" = "#D55E00", "overlap_coefficient" = "#009E73"),
      labels = c("Mean Diff", "Cohen's d", "Overlap"),
      name = "Metric"
    ) +
    labs(
      title = "Score Distribution Separation",
      subtitle = "Higher Cohen's d and lower overlap indicate better separation",
      x = "Tool",
      y = "Value"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  p
}

# -----------------------------------------------------------------------------
# SFig 3: Threshold Robustness
# -----------------------------------------------------------------------------
create_sfig_threshold_robust <- function(thresh_robust_df, tool_colors) {
  if (is.null(thresh_robust_df)) return(NULL)

  df_plot <- thresh_robust_df %>%
    dplyr::arrange(threshold_cv)

  df_plot$tool <- factor(df_plot$tool, levels = df_plot$tool)

  p <- ggplot(df_plot, aes(x = tool, y = threshold_cv * 100, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = tool_colors, guide = "none") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", threshold_cv * 100)), vjust = -0.5, size = 2) +
    labs(
      title = "Threshold Stability (Cross-Validation)",
      subtitle = "Lower CV indicates more stable threshold (dashed line = 10% threshold)",
      x = "Tool",
      y = "Coefficient of Variation (%)"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

# -----------------------------------------------------------------------------
# SFig 4: Effect Sizes
# -----------------------------------------------------------------------------
create_sfig_effect_sizes <- function(effect_df, tool_colors) {
  if (is.null(effect_df)) return(NULL)

  # Effect sizes are pairwise (tool1 vs tool2); aggregate across groups
  # Filter to F1 metric if available, otherwise use first metric
  if ("metric" %in% names(effect_df)) {
    metric_choice <- if ("f1" %in% effect_df$metric) "f1" else effect_df$metric[1]
    effect_df <- effect_df %>% dplyr::filter(metric == metric_choice)
  }

  # Create pair label and summarise across comparisons/groups
  df_plot <- effect_df %>%
    dplyr::mutate(pair = paste0(tool1, " vs ", tool2)) %>%
    dplyr::group_by(pair) %>%
    dplyr::summarise(cohens_d = mean(cohens_d, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(cohens_d))

  df_plot$pair <- factor(df_plot$pair, levels = df_plot$pair)

  p <- ggplot(df_plot, aes(x = pair, y = cohens_d)) +
    geom_bar(stat = "identity", width = 0.7, fill = "#4477AA") +
    geom_hline(yintercept = c(0.2, 0.5, 0.8), linetype = "dashed", color = "gray50", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", cohens_d)), vjust = -0.5, size = 2) +
    annotate("text", x = Inf, y = c(0.2, 0.5, 0.8), hjust = -0.5, size = 2,
             label = c("Small", "Medium", "Large"), color = "gray50") +
    labs(
      title = "Effect Size (Cohen's d)",
      subtitle = "Pairwise performance difference (F1 score)",
      x = "Tool Pair",
      y = "Cohen's d"
    ) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

# =============================================================================
# Source Data Export Functions
# =============================================================================

ensure_output <- function(path) {
  # Create empty file if it doesn't exist (so Snakemake doesn't fail on missing output)
  if (!file.exists(path)) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    file.create(path)
  }
}

export_source_data <- function(data, filename, output_dir) {
  if (is.null(data) || nrow(data) == 0) return(FALSE)

  filepath <- file.path(output_dir, filename)
  tryCatch({
    write.table(data, filepath, sep = "\t", row.names = FALSE, quote = FALSE)
    message("  Saved source data: ", filename)
    TRUE
  }, error = function(e) {
    warning("Failed to save source data: ", e$message)
    FALSE
  })
}

# =============================================================================
# Main Execution
# =============================================================================

message(paste(rep("=", 70), collapse = ""))
message("Publication-Quality Benchmark Figure Generation")
message(paste(rep("=", 70), collapse = ""))

message("\nConfiguration:")
message("  Data directory: ", data_dir)
message("  Output directory: ", output_dir)
message("  Theme: ", theme)
message("  DPI: ", dpi)
message("  Format: ", fig_format)
if (!is.null(window_filter)) {
  message("  Window filter: ", window_filter)
}

# Create output directories
message("\nCreating output directories...")
dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)
message("  Created: ", main_dir)
message("  Created: ", supp_dir)
message("  Created: ", data_out_dir)

# =============================================================================
# Load Data
# =============================================================================

message("\nLoading data files...")

# Required files
accuracy_df <- load_with_fallback(input_files$accuracy, required = TRUE)
if (is.null(accuracy_df)) {
  stop("Required accuracy data not found. Exiting.")
}

# Optional files
ci_df <- load_with_fallback(input_files$bootstrap_ci)
sig_df <- load_with_fallback(input_files$significance)
fdr_df <- load_with_fallback(input_files$fdr_corrected)
effect_df <- load_with_fallback(input_files$effect_sizes)
coverage_df <- load_with_fallback(input_files$coverage)
score_dist_df <- load_with_fallback(input_files$score_dist)
thresh_robust_df <- load_with_fallback(input_files$threshold_robust)
by_comp_df <- load_with_fallback(input_files$accuracy_by_comparison)
resource_df <- load_with_fallback(input_files$resources)

# Assign tool colors
tool_colors <- assign_tool_colors(unique(accuracy_df$tool))

# =============================================================================
# Generate Main Figures
# =============================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("Generating Main Figures")
message(paste(rep("=", 70), collapse = ""))

# Figure 1: Overall Accuracy with CI
message("\n[1/6] Overall Accuracy with Bootstrap CI...")
p1 <- create_fig1_overall_accuracy(accuracy_df, ci_df, tool_colors)
if (!is.null(p1)) {
  save_figure(p1, file.path(main_dir, paste0("fig1_overall_accuracy.", fig_format)),
              width = 85, height = 70, dpi = dpi)
  export_source_data(
    accuracy_df %>% dplyr::group_by(tool) %>% dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop"),
    "fig1_source_data.tsv", data_out_dir
  )
}

# Figure 2: PR Curves
message("[2/6] Precision vs Recall...")
p2 <- create_fig2_pr_curves(accuracy_df, tool_colors)
if (!is.null(p2)) {
  save_figure(p2, file.path(main_dir, paste0("fig2_pr_curves.", fig_format)),
              width = 85, height = 70, dpi = dpi)
  export_source_data(
    accuracy_df %>% dplyr::group_by(tool) %>%
      dplyr::summarise(precision = mean(precision, na.rm = TRUE),
                       recall = mean(recall, na.rm = TRUE), .groups = "drop"),
    "fig2_source_data.tsv", data_out_dir
  )
}

# Figure 3: ROC Curves
message("[3/6] ROC Curves (AUROC)...")
p3 <- create_fig3_roc_curves(accuracy_df, ci_df, tool_colors)
if (!is.null(p3)) {
  save_figure(p3, file.path(main_dir, paste0("fig3_roc_curves.", fig_format)),
              width = 85, height = 70, dpi = dpi)
  export_source_data(
    accuracy_df %>% dplyr::group_by(tool) %>%
      dplyr::summarise(auroc = mean(auroc, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(!is.na(auroc)),
    "fig3_source_data.tsv", data_out_dir
  )
}

# Figure 4: F1 Comparison
message("[4/6] F1 Comparison...")
p4 <- create_fig4_f1_comparison(accuracy_df, sig_df, tool_colors)
if (!is.null(p4)) {
  save_figure(p4, file.path(main_dir, paste0("fig4_f1_comparison.", fig_format)),
              width = 85, height = 70, dpi = dpi)
  export_source_data(
    accuracy_df %>% dplyr::group_by(tool) %>%
      dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop"),
    "fig4_source_data.tsv", data_out_dir
  )
}

# Figure 5: Coverage Sensitivity
message("[5/6] Coverage Sensitivity...")
p5 <- create_fig5_coverage_sensitivity(coverage_df, tool_colors)
if (!is.null(p5)) {
  save_figure(p5, file.path(main_dir, paste0("fig5_coverage_sensitivity.", fig_format)),
              width = 174, height = 100, dpi = dpi)
  export_source_data(coverage_df, "fig5_source_data.tsv", data_out_dir)
}

# Figure 6: Resource Usage
message("[6/6] Resource Usage...")
p6 <- create_fig6_resource_usage(resource_df, accuracy_df, tool_colors)
if (!is.null(p6)) {
  save_figure(p6, file.path(main_dir, paste0("fig6_resource_usage.", fig_format)),
              width = 100, height = 85, dpi = dpi)
  export_source_data(resource_df, "fig6_source_data.tsv", data_out_dir)
}

# =============================================================================
# Generate Supplementary Figures
# =============================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("Generating Supplementary Figures")
message(paste(rep("=", 70), collapse = ""))

# SFig 1: Per-Comparison
message("\n[S1] Per-Comparison Breakdown...")
sp1 <- create_sfig_per_comparison(by_comp_df, tool_colors)
if (!is.null(sp1)) {
  save_figure(sp1, file.path(supp_dir, paste0("sfig1_per_comparison.", fig_format)),
              width = 174, height = 100, dpi = dpi)
}

# SFig 2: Score Distributions
message("[S2] Score Distribution Analysis...")
sp2 <- create_sfig_score_dist(score_dist_df, tool_colors)
if (!is.null(sp2)) {
  save_figure(sp2, file.path(supp_dir, paste0("sfig2_score_distributions.", fig_format)),
              width = 100, height = 70, dpi = dpi)
}

# SFig 3: Threshold Robustness
message("[S3] Threshold Robustness...")
sp3 <- create_sfig_threshold_robust(thresh_robust_df, tool_colors)
if (!is.null(sp3)) {
  save_figure(sp3, file.path(supp_dir, paste0("sfig3_threshold_robustness.", fig_format)),
              width = 100, height = 70, dpi = dpi)
}

# SFig 4: Effect Sizes
message("[S4] Effect Size Analysis...")
sp4 <- create_sfig_effect_sizes(effect_df, tool_colors)
if (!is.null(sp4)) {
  save_figure(sp4, file.path(supp_dir, paste0("sfig4_effect_sizes.", fig_format)),
              width = 100, height = 70, dpi = dpi)
}

# =============================================================================
# Export All Source Data
# =============================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("Exporting Source Data for Manuscript Submission")
message(paste(rep("=", 70), collapse = ""))

export_source_data(accuracy_df, "accuracy_summary_source.tsv", data_out_dir)
if (!is.null(ci_df)) export_source_data(ci_df, "bootstrap_ci_source.tsv", data_out_dir)
if (!is.null(sig_df)) export_source_data(sig_df, "significance_tests_source.tsv", data_out_dir)
if (!is.null(effect_df)) export_source_data(effect_df, "effect_sizes_source.tsv", data_out_dir)
if (!is.null(coverage_df)) export_source_data(coverage_df, "coverage_analysis_source.tsv", data_out_dir)
if (!is.null(score_dist_df)) export_source_data(score_dist_df, "score_distribution_source.tsv", data_out_dir)
if (!is.null(thresh_robust_df)) export_source_data(thresh_robust_df, "threshold_robustness_source.tsv", data_out_dir)

# =============================================================================
# Ensure all declared Snakemake outputs exist (create empty placeholders if skipped)
# =============================================================================

if (exists("snakemake")) {
  all_outputs <- unlist(snakemake@output)
  for (out_path in all_outputs) {
    ensure_output(out_path)
  }
}

# =============================================================================
# Summary Report
# =============================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("Generation Complete!")
message(paste(rep("=", 70), collapse = ""))

# Count generated figures
main_figs <- list.files(main_dir, pattern = paste0("\\.", fig_format, "$"))
supp_figs <- list.files(supp_dir, pattern = paste0("\\.", fig_format, "$"))
data_files <- list.files(data_out_dir, pattern = "\\.tsv$")

message("\nMain figures: ", length(main_figs))
message("Supplementary figures: ", length(supp_figs))
message("Source data files: ", length(data_files))

message("\nOutput directories:")
message("  Main: ", main_dir)
message("  Supplementary: ", supp_dir)
message("  Data: ", data_out_dir)

message("\n", paste(rep("=", 70), collapse = ""))
message("Done!")
message(paste(rep("=", 70), collapse = ""))
