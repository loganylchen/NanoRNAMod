#!/usr/bin/env Rscript
# =============================================================================
# Per-Modification Analysis Figure
# =============================================================================
# Generates stratified visualizations by modification type:
# - Heatmap of F1 scores (tools × modification types)
# - Grouped bar chart showing performance per modification type
# - Faceted plots for each major modification type
#
# Usage:
#   Rscript 02_per_modification_analysis.R --input accuracy_summary.tsv --output figures/
# =============================================================================

options(echo = FALSE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- "accuracy_summary.tsv"
output_dir <- "figures"

for (i in seq_along(args)) {
  if (args[i] == "--input" && i + 1 <= length(args)) {
    input_file <- args[i + 1]
  } else if (args[i] == "--output" && i + 1 <= length(args)) {
    output_dir <- args[i + 1]
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
# Data Loading and Preparation
# =============================================================================

message("Loading data from: ", input_file)
df <- load_accuracy_summary(input_file)

# Check if modification_type column exists
if (!"modification_type" %in% names(df)) {
  stop("modification_type column not found in data. ",
       "This script requires per-modification data.")
}

# Filter out "overall" category for stratified analysis
df_mod <- df[df$modification_type != "overall" & !is.na(df$modification_type), ]

if (nrow(df_mod) == 0) {
  stop("No per-modification data found after filtering.")
}

message("Found modification types: ", paste(unique(df_mod$modification_type), collapse = ", "))

# Assign colors
mod_colors <- okabe_ito_palette(length(unique(df_mod$modification_type)))
names(mod_colors) <- unique(df_mod$modification_type)

# =============================================================================
# Figure 1: F1 Score Heatmap (Tools × Modifications)
# =============================================================================

create_f1_heatmap <- function(data) {
  # Aggregate by tool and modification
  df_heat <- data %>%
    dplyr::group_by(tool, modification_type) %>%
    dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(tool, modification_type, fill = list(f1 = NA_real_)) %>%
    tidyr::pivot_wider(names_from = modification_type, values_from = f1)

  # Set row names for pheatmap
  df_matrix <- as.matrix(df_heat[, -1])
  rownames(df_matrix) <- df_heat$tool

  # Reorder by mean F1
  row_order <- rowMeans(df_matrix, na.rm = TRUE) %>%
    order(decreasing = TRUE)
  df_matrix <- df_matrix[row_order, ]

  # Column order by mean F1
  col_order <- colMeans(df_matrix, na.rm = TRUE) %>%
    order(decreasing = TRUE)
  df_matrix <- df_matrix[, col_order]

  # Create heatmap using ggplot2
  df_long <- data %>%
    dplyr::group_by(tool, modification_type) %>%
    dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop")

  df_long$tool <- factor(df_long$tool, levels = rownames(df_matrix))
  df_long$modification_type <- factor(df_long$modification_type,
                                      levels = colnames(df_matrix))

  p <- ggplot(df_long, aes(x = modification_type, y = tool, fill = f1)) +
    geom_tile(color = "white", size = 0.25) +
    geom_text(aes(label = ifelse(!is.na(f1), sprintf("%.2f", f1), "")),
              size = 2, color = "black") +
    scale_fill_gradient2(
      low = "#D73027", mid = "#FFFFBF", high = "#1A9850",
      midpoint = 0.5, na.value = "gray90",
      limits = c(0, 1), name = "F1 Score"
    ) +
    labs(
      title = "F1 Score by Tool and Modification Type",
      x = "Modification Type",
      y = "Tool"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      legend.position = "right"
    )

  p
}

# =============================================================================
# Figure 2: Grouped Bar Chart (Per Modification Type)
# =============================================================================

create_grouped_bar <- function(data, colors) {
  df_summary <- data %>%
    dplyr::group_by(tool, modification_type) %>%
    dplyr::summarise(
      f1 = mean(f1, na.rm = TRUE),
      precision = mean(precision, na.rm = TRUE),
      recall = mean(recall, na.rm = TRUE),
      .groups = "drop"
    )

  # Order tools by overall F1
  tool_order <- data %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(f1)) %>%
    pull(tool)

  df_summary$tool <- factor(df_summary$tool, levels = tool_order)

  p <- ggplot(df_summary, aes(x = tool, y = f1, fill = modification_type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9),
             width = 0.8, color = "white", size = 0.25) +
    scale_fill_manual(values = colors) +
    labs(
      title = "F1 Score by Tool and Modification Type",
      subtitle = "Grouped comparison across modification types",
      x = "Tool",
      y = "F1 Score",
      fill = "Modification"
    ) +
    scale_y_continuous(limits = c(0, min(1.1, max(df_summary$f1, na.rm = TRUE) * 1.1)),
                       expand = c(0, 0)) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  p
}

# =============================================================================
# Figure 3: Faceted Precision-Recall by Modification
# =============================================================================

create_faceted_pr <- function(data) {
  df_summary <- data %>%
    dplyr::group_by(tool, modification_type) %>%
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

  # Order tools
  tool_order <- data %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(f1)) %>%
    pull(tool)

  df_summary$tool <- factor(df_summary$tool, levels = tool_order)

  p <- ggplot(df_summary, aes(x = tool, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7) +
    scale_fill_manual(
      values = c("precision" = "#0072B2", "recall" = "#D55E00"),
      labels = c("Precision", "Recall")
    ) +
    facet_wrap(~ modification_type, ncol = 2) +
    labs(
      title = "Precision and Recall by Modification Type",
      x = "Tool",
      y = "Score",
      fill = "Metric"
    ) +
    scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      strip.text = element_text(size = 7),
      legend.position = "top"
    )

  p
}

# =============================================================================
# Figure 4: Modification-Specific Faceted F1 Comparison
# =============================================================================

create_modification_facets <- function(data) {
  # Calculate summary for each modification type
  df_summary <- data %>%
    dplyr::group_by(tool, modification_type, window) %>%
    dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop")

  # Overall tool ordering
  tool_order <- data %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(f1)) %>%
    pull(tool)

  df_summary$tool <- factor(df_summary$tool, levels = tool_order)

  p <- ggplot(df_summary, aes(x = window, y = f1, color = tool)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.5) +
    scale_color_manual(values = assign_tool_colors(tool_order)) +
    facet_wrap(~ modification_type, ncol = 2) +
    labs(
      title = "F1 Score Across Windows by Modification Type",
      subtitle = "Each panel shows a different modification type",
      x = "Window (nt)",
      y = "F1 Score",
      color = "Tool"
    ) +
    scale_y_continuous(limits = c(0, 1.0), expand = c(0, 0)) +
    theme_nature() +
    theme(
      strip.text = element_text(size = 7),
      legend.position = "right",
      axis.text = element_text(size = 6),
      legend.text = element_text(size = 6)
    )

  p
}

# =============================================================================
# Figure 5: Stacked Bar of Modification Coverage
# =============================================================================

create_coverage_bar <- function(data) {
  # Count predictions per modification
  df_counts <- data %>%
    dplyr::group_by(tool, modification_type) %>%
    dplyr::summarise(
      tp = sum(tp, na.rm = TRUE),
      fp = sum(fp, na.rm = TRUE),
      fn = sum(fn, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(cols = c(tp, fp, fn), names_to = "count_type",
                       values_to = "count")

  # Order by total calls
  tool_order <- df_counts %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(total = sum(count, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(total)) %>%
    pull(tool)

  df_counts$tool <- factor(df_counts$tool, levels = tool_order)

  p <- ggplot(df_counts, aes(x = tool, y = count, fill = count_type)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7,
             color = "white", size = 0.25) +
    scale_fill_manual(
      values = c("tp" = "#009E73", "fp" = "#E69F00", "fn" = "#D55E00"),
      labels = c("TP", "FP", "FN"),
      name = "Count Type"
    ) +
    facet_wrap(~ modification_type, scales = "free_y") +
    labs(
      title = "Prediction Counts by Tool and Modification Type",
      x = "Tool",
      y = "Count"
    ) +
    scale_y_continuous(labels = scales::comma_format()) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      strip.text = element_text(size = 7),
      legend.position = "top"
    )

  p
}

# =============================================================================
# Generate and Save Figures
# =============================================================================

message("Generating per-modification figures...")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Figure 1: Heatmap
p1 <- create_f1_heatmap(df_mod)
save_figure(p1, file.path(output_dir, "modification_01_heatmap.pdf"),
            width = 120, height = 100)
save_figure(p1, file.path(output_dir, "modification_01_heatmap.png"),
            width = 120, height = 100, dpi = 300)

# Figure 2: Grouped bar
p2 <- create_grouped_bar(df_mod, mod_colors)
save_figure(p2, file.path(output_dir, "modification_02_grouped_bar.pdf"),
            width = 150, height = 85)
save_figure(p2, file.path(output_dir, "modification_02_grouped_bar.png"),
            width = 150, height = 85, dpi = 300)

# Figure 3: Faceted PR
p3 <- create_faceted_pr(df_mod)
save_figure(p3, file.path(output_dir, "modification_03_faceted_pr.pdf"),
            width = 120, height = 120)
save_figure(p3, file.path(output_dir, "modification_03_faceted_pr.png"),
            width = 120, height = 120, dpi = 300)

# Figure 4: Window faceted
if ("window" %in% names(df_mod) && length(unique(df_mod$window)) > 1) {
  p4 <- create_modification_facets(df_mod)
  save_figure(p4, file.path(output_dir, "modification_04_window_facets.pdf"),
              width = 140, height = 120)
  save_figure(p4, file.path(output_dir, "modification_04_window_facets.png"),
              width = 140, height = 120, dpi = 300)
}

# Figure 5: Coverage
if (all(c("tp", "fp", "fn") %in% names(df_mod))) {
  p5 <- create_coverage_bar(df_mod)
  save_figure(p5, file.path(output_dir, "modification_05_coverage.pdf"),
              width = 140, height = 100)
  save_figure(p5, file.path(output_dir, "modification_05_coverage.png"),
              width = 140, height = 100, dpi = 300)
}

message("Done! Per-modification figures saved to: ", output_dir)
