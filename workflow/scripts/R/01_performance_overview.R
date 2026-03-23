#!/usr/bin/env Rscript
# =============================================================================
# Performance Overview Figure
# =============================================================================
# Generates main performance comparison visualizations:
# - F1 score bar chart (all tools sorted by performance)
# - Precision and Recall grouped bar chart
# - ROC curves (if AUROC data available)
# - PR curves (if AUPRC data available)
#
# Usage:
#   Rscript 01_performance_overview.R --input accuracy_summary.tsv --output figures/
# =============================================================================

# Suppress startup messages
options(echo = FALSE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values
input_file <- "accuracy_summary.tsv"
output_dir <- "figures"
window_filter <- NULL

# Parse arguments
if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--input" && i + 1 <= length(args)) {
      input_file <- args[i + 1]
    } else if (args[i] == "--output" && i + 1 <= length(args)) {
      output_dir <- args[i + 1]
    } else if (args[i] == "--window" && i + 1 <= length(args)) {
      window_filter <- as.integer(args[i + 1])
    }
  }
}

# Load utilities
script_dir <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      path <- sub("^--file=", "", file_arg[1])
      dirname(normalizePath(path))
    } else {
      dirname(sys.frame(1)$ofile)
    }
  }, error = function(e) ".")
source(file.path(script_dir, "00_utils.R"))

# =============================================================================
# Main Analysis
# =============================================================================

message("Loading data from: ", input_file)
df <- load_accuracy_summary(input_file)

# Filter by window if specified
if (!is.null(window_filter)) {
  df <- df[df$window == window_filter, ]
  message("Filtered to window: ", window_filter)
}

# Filter to overall summary (no modification_type) if available
if ("modification_type" %in% names(df)) {
  df_overall <- df[df$modification_type == "overall" | is.na(df$modification_type), ]
  if (nrow(df_overall) > 0) {
    df <- df_overall
    message("Using overall summary (aggregated across modifications)")
  }
}

# Assign tool colors
tool_colors <- assign_tool_colors(unique(df$tool))

# =============================================================================
# Figure 1: F1 Score Bar Chart (Sorted)
# =============================================================================

create_f1_bar_chart <- function(data, colors) {
  # Aggregate by tool (mean across windows/modifications)
  df_summary <- data %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(
      f1 = mean(f1, na.rm = TRUE),
      precision = mean(precision, na.rm = TRUE),
      recall = mean(recall, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(f1))

  # Reorder factor levels
  df_summary$tool <- factor(df_summary$tool, levels = df_summary$tool)

  # Match colors to tool order
  plot_colors <- colors[match(levels(df_summary$tool), names(colors))]

  p <- ggplot(df_summary, aes(x = tool, y = f1, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = plot_colors, guide = "none") +
    geom_text(aes(label = sprintf("%.2f", f1)),
              vjust = -0.5, size = 2) +
    labs(
      title = "Tool Performance Comparison",
      subtitle = "F1 Score (higher is better)",
      x = "Tool",
      y = "F1 Score"
    ) +
    scale_y_continuous(limits = c(0, max(1.1, max(df_summary$f1) * 1.1)),
                       expand = c(0, 0)) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

# =============================================================================
# Figure 2: Precision and Recall Grouped Bar Chart
# =============================================================================

create_precision_recall_bar <- function(data, colors) {
  df_long <- data %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(
      precision = mean(precision, na.rm = TRUE),
      recall = mean(recall, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(precision, recall),
      names_to = "metric",
      values_to = "value"
    )

  # Order by F1 score for tool arrangement
  tool_order <- data %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(f1)) %>%
    pull(tool)

  df_long$tool <- factor(df_long$tool, levels = tool_order)

  p <- ggplot(df_long, aes(x = tool, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7) +
    scale_fill_manual(
      values = c("precision" = "#0072B2", "recall" = "#D55E00"),
      labels = c("Precision", "Recall")
    ) +
    geom_text(aes(label = sprintf("%.2f", value)),
              position = position_dodge(width = 0.8),
              vjust = -0.5, size = 2) +
    labs(
      title = "Precision vs Recall",
      subtitle = "Trade-off between precision and recall",
      x = "Tool",
      y = "Score"
    ) +
    scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")

  p
}

# =============================================================================
# Figure 3: ROC Curves (AUROC Bar Chart)
# =============================================================================

create_auroc_bar <- function(data, colors) {
  if (!"auroc" %in% names(data) || all(is.na(data$auroc))) {
    message("AUROC data not available")
    return(NULL)
  }

  df_summary <- data %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(auroc = mean(auroc, na.rm = TRUE), .groups = "drop") %>%
    tidyr::drop_na() %>%
    dplyr::arrange(desc(auroc))

  df_summary$tool <- factor(df_summary$tool, levels = df_summary$tool)
  plot_colors <- colors[match(levels(df_summary$tool), names(colors))]

  p <- ggplot(df_summary, aes(x = tool, y = auroc, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = plot_colors, guide = "none") +
    geom_text(aes(label = sprintf("%.3f", auroc)),
              vjust = -0.5, size = 2) +
    labs(
      title = "AUROC Comparison",
      subtitle = "Area Under ROC Curve (higher is better)",
      x = "Tool",
      y = "AUROC"
    ) +
    scale_y_continuous(limits = c(0.5, 1.0), expand = c(0, 0)) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "gray50", linewidth = 0.25) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

# =============================================================================
# Figure 4: PR Curves (AUPRC Bar Chart)
# =============================================================================

create_auprc_bar <- function(data, colors) {
  if (!"auprc" %in% names(data) || all(is.na(data$auprc))) {
    message("AUPRC data not available")
    return(NULL)
  }

  df_summary <- data %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(auprc = mean(auprc, na.rm = TRUE), .groups = "drop") %>%
    tidyr::drop_na() %>%
    dplyr::arrange(desc(auprc))

  df_summary$tool <- factor(df_summary$tool, levels = df_summary$tool)
  plot_colors <- colors[match(levels(df_summary$tool), names(colors))]

  p <- ggplot(df_summary, aes(x = tool, y = auprc, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = plot_colors, guide = "none") +
    geom_text(aes(label = sprintf("%.3f", auprc)),
              vjust = -0.5, size = 2) +
    labs(
      title = "AUPRC Comparison",
      subtitle = "Area Under PR Curve (higher is better)",
      x = "Tool",
      y = "AUPRC"
    ) +
    scale_y_continuous(limits = c(0, max(1.1, max(df_summary$auprc) * 1.1)),
                       expand = c(0, 0)) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

# =============================================================================
# Figure 5: Combined Metrics Radar/Spider Plot
# =============================================================================

create_metrics_spider <- function(data) {
  metrics <- c("f1", "precision", "recall")
  available <- intersect(metrics, names(data))

  if (length(available) < 2) {
    message("Not enough metrics for spider plot")
    return(NULL)
  }

  df_summary <- data %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(across(all_of(available), ~mean(.x, na.rm = TRUE)),
                    .groups = "drop")

  # Skip if ggradar not available
  if (!requireNamespace("ggradar", quietly = TRUE)) {
    message("ggradar package not available for spider plot")
    return(NULL)
  }

  # Normalize to 0-1 scale for radar
  for (col in available) {
    max_val <- max(df_summary[[col]], na.rm = TRUE)
    if (max_val > 0) {
      df_summary[[col]] <- df_summary[[col]] / max_val
    }
  }

  p <- ggradar::ggradar(
    df_summary,
    values.radar = c(0, 1),
    grid.min = 0, grid.mid = 0.5, grid.max = 1,
    colours = okabe_ito_palette()
  ) +
    labs(title = "Multi-Metric Comparison")

  p
}

# =============================================================================
# Generate and Save Figures
# =============================================================================

message("Generating figures...")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Figure 1: F1 Score
p1 <- create_f1_bar_chart(df, tool_colors)
save_figure(p1, file.path(output_dir, "01_f1_scores.pdf"),
            width = 85, height = 70)
save_figure(p1, file.path(output_dir, "01_f1_scores.png"),
            width = 85, height = 70, dpi = 300)

# Figure 2: Precision vs Recall
p2 <- create_precision_recall_bar(df, tool_colors)
save_figure(p2, file.path(output_dir, "02_precision_recall.pdf"),
            width = 85, height = 70)
save_figure(p2, file.path(output_dir, "02_precision_recall.png"),
            width = 85, height = 70, dpi = 300)

# Figure 3: AUROC (if available)
p3 <- create_auroc_bar(df, tool_colors)
if (!is.null(p3)) {
  save_figure(p3, file.path(output_dir, "03_auroc.pdf"),
              width = 85, height = 70)
  save_figure(p3, file.path(output_dir, "03_auroc.png"),
              width = 85, height = 70, dpi = 300)
}

# Figure 4: AUPRC (if available)
p4 <- create_auprc_bar(df, tool_colors)
if (!is.null(p4)) {
  save_figure(p4, file.path(output_dir, "04_auprc.pdf"),
              width = 85, height = 70)
  save_figure(p4, file.path(output_dir, "04_auprc.png"),
              width = 85, height = 70, dpi = 300)
}

# Figure 5: Spider plot (optional)
p5 <- create_metrics_spider(df)
if (!is.null(p5)) {
  save_figure(p5, file.path(output_dir, "05_metrics_spider.pdf"),
              width = 100, height = 100)
}

message("Done! Figures saved to: ", output_dir)
