#!/usr/bin/env Rscript
# =============================================================================
# Threshold Optimization Analysis Figure
# =============================================================================
# Generates threshold analysis visualizations:
# - Line plots showing precision-recall tradeoff at different thresholds
# - Optimal threshold annotation on curves
# - Faceted by tool
#
# Usage:
#   Rscript 03_threshold_optimization.R --input threshold_evaluation.tsv --output figures/
# =============================================================================

options(echo = FALSE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- "threshold_evaluation.tsv"
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
# Data Loading
# =============================================================================

message("Loading threshold data from: ", input_file)

if (!file.exists(input_file)) {
  stop("Threshold evaluation file not found: ", input_file)
}

df <- readr::read_tsv(input_file, show_col_types = FALSE)

message("Data dimensions: ", nrow(df), " rows, ", ncol(df), " columns")
message("Columns: ", paste(names(df), collapse = ", "))

# Check required columns
required_cols <- c("tool", "threshold", "precision", "recall", "f1")
missing <- setdiff(required_cols, names(df))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

# Assign tool colors
tool_colors <- assign_tool_colors(unique(df$tool))

# =============================================================================
# Figure 1: Precision-Recall Tradeoff by Tool
# =============================================================================

create_pr_tradeoff <- function(data, colors) {
  # Aggregate thresholds (some tools may have multiple score columns)
  df_summary <- data %>%
    dplyr::group_by(tool, threshold) %>%
    dplyr::summarise(
      precision = mean(precision, na.rm = TRUE),
      recall = mean(recall, na.rm = TRUE),
      f1 = mean(f1, na.rm = TRUE),
      .groups = "drop"
    )

  # Find optimal threshold per tool (max F1)
  optimal <- df_summary %>%
    dplyr::group_by(tool) %>%
    dplyr::slice(which.max(f1)) %>%
    dplyr::ungroup()

  p <- ggplot(df_summary, aes(x = threshold, y = precision, color = tool)) +
    geom_line(linewidth = 0.4, alpha = 0.8) +
    geom_point(data = optimal,
               aes(x = threshold, y = precision),
               size = 2, shape = 16) +
    scale_color_manual(values = colors) +
    labs(
      title = "Precision-Recall Tradeoff by Threshold",
      subtitle = "Points indicate optimal threshold (max F1)",
      x = "Threshold",
      y = "Precision",
      color = "Tool"
    ) +
    theme_nature() +
    theme(legend.position = "right")

  p
}

# =============================================================================
# Figure 2: Precision and Recall Curves (Faceted)
# =============================================================================

create_pr_curves_faceted <- function(data) {
  df_long <- data %>%
    dplyr::group_by(tool, threshold, score_column) %>%
    dplyr::summarise(
      precision = mean(precision, na.rm = TRUE),
      recall = mean(recall, na.rm = TRUE),
      f1 = mean(f1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(precision, recall, f1),
      names_to = "metric",
      values_to = "value"
    )

  # Find optimal points
  optimal <- df_long %>%
    dplyr::filter(metric == "f1") %>%
    dplyr::group_by(tool, score_column) %>%
    dplyr::slice(which.max(value)) %>%
    dplyr::ungroup() %>%
    dplyr::select(tool, score_column, threshold, optimal_f1 = value)

  # Merge optimal thresholds back
  df_long <- df_long %>%
    dplyr::left_join(optimal,
                     by = c("tool", "score_column", "threshold"),
                     suffix = c("", "_opt"))

  # Order tools by overall performance
  tool_order <- data %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(f1 = max(f1, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(f1)) %>%
    pull(tool)

  df_long$tool <- factor(df_long$tool, levels = tool_order)

  p <- ggplot(df_long, aes(x = threshold, y = value, color = metric)) +
    geom_line(linewidth = 0.5) +
    geom_point(data = df_long[!is.na(df_long$optimal_f1), ],
               shape = 16, size = 1.5) +
    scale_color_manual(
      values = c("precision" = "#0072B2",
                 "recall" = "#D55E00",
                 "f1" = "#009E73"),
      labels = c("Precision", "Recall", "F1")
    ) +
    facet_wrap(~ tool, ncol = 3) +
    labs(
      title = "Threshold Sweep: Precision, Recall, and F1",
      subtitle = "Dots indicate optimal threshold for each tool",
      x = "Threshold",
      y = "Score",
      color = "Metric"
    ) +
    theme_nature() +
    theme(
      strip.text = element_text(size = 7),
      legend.position = "top",
      axis.text = element_text(size = 6)
    )

  p
}

# =============================================================================
# Figure 3: F1 vs Threshold with Optimal Highlighted
# =============================================================================

create_f1_threshold_lines <- function(data, colors) {
  # Handle multiple score columns per tool
  df_summary <- data %>%
    dplyr::group_by(tool, score_column, threshold) %>%
    dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop")

  # Create unique identifier
  df_summary$tool_score <- ifelse(
    is.na(df_summary$score_column) | df_summary$score_column == "",
    df_summary$tool,
    paste(df_summary$tool, df_summary$score_column, sep = ": ")
  )

  # Find optimal
  optimal <- df_summary %>%
    dplyr::group_by(tool, score_column) %>%
    dplyr::slice(which.max(f1)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(tool_score = ifelse(
      is.na(score_column) | score_column == "",
      tool,
      paste(tool, score_column, sep = ": ")
    ))

  # Order by max F1
  item_order <- df_summary %>%
    dplyr::group_by(tool, score_column) %>%
    dplyr::summarise(max_f1 = max(f1, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(max_f1)) %>%
    dplyr::mutate(tool_score = ifelse(
      is.na(score_column) | score_column == "",
      tool,
      paste(tool, score_column, sep = ": ")
    )) %>%
    pull(tool_score)

  df_summary$tool_score <- factor(df_summary$tool_score, levels = item_order)

  # Assign colors
  n_items <- length(unique(df_summary$tool_score))
  item_colors <- rep(okabe_ito_palette(), length.out = n_items)
  names(item_colors) <- levels(df_summary$tool_score)

  p <- ggplot(df_summary, aes(x = threshold, y = f1, color = tool_score)) +
    geom_line(linewidth = 0.5) +
    geom_point(data = optimal, shape = 16, size = 2) +
    scale_color_manual(values = item_colors) +
    labs(
      title = "F1 Score vs Threshold",
      subtitle = "Each line represents a tool/score column combination. Dots = optimal.",
      x = "Threshold",
      y = "F1 Score",
      color = "Tool"
    ) +
    theme_nature() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 6)
    )

  p
}

# =============================================================================
# Figure 4: Threshold Stability Analysis
# =============================================================================

create_threshold_stability <- function(data) {
  # Calculate coefficient of variation across thresholds near optimum
  df_summary <- data %>%
    dplyr::group_by(tool, score_column, threshold) %>%
    dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop")

  # Define "near optimal" as within 10% of optimal threshold range
  stability <- df_summary %>%
    dplyr::group_by(tool, score_column) %>%
    dplyr::mutate(
      max_f1 = max(f1, na.rm = TRUE),
      is_near_optimal = f1 >= 0.9 * max_f1
    ) %>%
    dplyr::filter(is_near_optimal) %>%
    dplyr::summarise(
      mean_f1_near_opt = mean(f1, na.rm = TRUE),
      sd_f1_near_opt = sd(f1, na.rm = TRUE),
      cv_f1 = ifelse(mean_f1_near_opt > 0, sd_f1_near_opt / mean_f1_near_opt, NA),
      n_thresholds_near_opt = dplyr::n(),
      optimal_f1 = max(f1, na.rm = TRUE),
      .groups = "drop"
    )

  # Create display name
  stability$tool_score <- ifelse(
    is.na(stability$score_column) | stability$score_column == "",
    stability$tool,
    paste(stability$tool, stability$score_column, sep = ": ")
  )

  stability <- stability %>%
    dplyr::arrange(desc(optimal_f1))

  stability$tool_score <- factor(stability$tool_score,
                                  levels = stability$tool_score)

  p <- ggplot(stability,
              aes(x = tool_score, y = cv_f1, fill = optimal_f1)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_gradient2(
      low = "#D73027", mid = "#FFFFBF", high = "#1A9850",
      midpoint = 0.7, name = "Optimal F1"
    ) +
    coord_flip() +
    labs(
      title = "Threshold Stability Analysis",
      subtitle = "Lower CV = more stable performance near optimal threshold",
      x = "Tool",
      y = "Coefficient of Variation\n(near optimal threshold)"
    ) +
    theme_nature()

  p
}

# =============================================================================
# Figure 5: Optimal Threshold Distribution
# =============================================================================

create_threshold_distribution <- function(data) {
  optimal <- data %>%
    dplyr::group_by(tool, score_column, threshold) %>%
    dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(tool, score_column) %>%
    dplyr::slice(which.max(f1)) %>%
    dplyr::ungroup()

  optimal$tool_score <- ifelse(
    is.na(optimal$score_column) | optimal$score_column == "",
    optimal$tool,
    paste(optimal$tool, optimal$score_column, sep = ": ")
  )

  # Order by threshold value
  optimal <- optimal %>%
    dplyr::arrange(threshold)

  optimal$tool_score <- factor(optimal$tool_score,
                                levels = optimal$tool_score)

  p <- ggplot(optimal, aes(x = tool_score, y = threshold, fill = f1)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_gradient(
      low = "#E69F00", high = "#009E73",
      name = "F1 Score"
    ) +
    coord_flip() +
    labs(
      title = "Optimal Threshold Distribution",
      subtitle = "Each bar shows the threshold that maximizes F1 for each tool",
      x = "Tool",
      y = "Optimal Threshold"
    ) +
    theme_nature()

  p
}

# =============================================================================
# Generate and Save Figures
# =============================================================================

message("Generating threshold optimization figures...")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Figure 1: PR Tradeoff
p1 <- create_pr_tradeoff(df, tool_colors)
save_figure(p1, file.path(output_dir, "threshold_01_pr_tradeoff.pdf"),
            width = 120, height = 85)
save_figure(p1, file.path(output_dir, "threshold_01_pr_tradeoff.png"),
            width = 120, height = 85, dpi = 300)

# Figure 2: Faceted PR curves
p2 <- create_pr_curves_faceted(df)
save_figure(p2, file.path(output_dir, "threshold_02_faceted_curves.pdf"),
            width = 174, height = 120)
save_figure(p2, file.path(output_dir, "threshold_02_faceted_curves.png"),
            width = 174, height = 120, dpi = 300)

# Figure 3: F1 lines
p3 <- create_f1_threshold_lines(df, tool_colors)
save_figure(p3, file.path(output_dir, "threshold_03_f1_lines.pdf"),
            width = 120, height = 85)
save_figure(p3, file.path(output_dir, "threshold_03_f1_lines.png"),
            width = 120, height = 85, dpi = 300)

# Figure 4: Stability
p4 <- create_threshold_stability(df)
save_figure(p4, file.path(output_dir, "threshold_04_stability.pdf"),
            width = 85, height = 100)
save_figure(p4, file.path(output_dir, "threshold_04_stability.png"),
            width = 85, height = 100, dpi = 300)

# Figure 5: Threshold distribution
p5 <- create_threshold_distribution(df)
save_figure(p5, file.path(output_dir, "threshold_05_distribution.pdf"),
            width = 85, height = 100)
save_figure(p5, file.path(output_dir, "threshold_05_distribution.png"),
            width = 85, height = 100, dpi = 300)

message("Done! Threshold optimization figures saved to: ", output_dir)
