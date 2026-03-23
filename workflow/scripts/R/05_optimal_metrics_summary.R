#!/usr/bin/env Rscript
# =============================================================================
# Optimal Metrics Summary Figure
# =============================================================================
# Generates summary visualization showing optimal configurations:
# - Table grob showing optimal score column and threshold per tool
# - Visual summary of which metric each tool should use
# - Summary statistics table
#
# Usage:
#   Rscript 05_optimal_metrics_summary.R --optimal optimal_thresholds.tsv --accuracy accuracy_summary.tsv --output figures/
# =============================================================================

options(echo = FALSE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
optimal_file <- "optimal_thresholds.tsv"
accuracy_file <- "accuracy_summary.tsv"
output_dir <- "figures"

for (i in seq_along(args)) {
  if (args[i] == "--optimal" && i + 1 <= length(args)) {
    optimal_file <- args[i + 1]
  } else if (args[i] == "--accuracy" && i + 1 <= length(args)) {
    accuracy_file <- args[i + 1]
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

message("Loading optimal thresholds from: ", optimal_file)

if (!file.exists(optimal_file)) {
  message("Optimal thresholds file not found. Will create summary from accuracy data.")
  optimal_df <- NULL
} else {
  optimal_df <- readr::read_tsv(optimal_file, show_col_types = FALSE)
}

message("Loading accuracy data from: ", accuracy_file)
accuracy_df <- load_accuracy_summary(accuracy_file)

# =============================================================================
# Summary Table Creation
# =============================================================================

create_optimal_summary_table <- function(optimal_data, accuracy_data) {
  # Build summary from optimal data
  if (!is.null(optimal_data) && nrow(optimal_data) > 0) {
    summary_df <- optimal_data %>%
      dplyr::select(tool, score_column, threshold,
                   any_of(c("f1", "precision", "recall"))) %>%
      dplyr::arrange(desc(f1))

    # Format for display
    summary_df <- summary_df %>%
      dplyr::mutate(
        score_label = ifelse(
          is.na(score_column) | score_column == "",
          "default",
          score_column
        ),
        threshold_fmt = sprintf("%.3f", threshold)
      )

    # Create display data frame
    display_df <- data.frame(
      Tool = summary_df$tool,
      `Score Column` = summary_df$score_label,
      `Optimal Threshold` = summary_df$threshold_fmt,
      F1 = sprintf("%.3f", summary_df$f1),
      Precision = sprintf("%.3f", summary_df$precision),
      Recall = sprintf("%.3f", summary_df$recall),
      check.names = FALSE
    )

    display_df
  } else {
    # Create summary from accuracy data (best window per tool)
    summary_df <- accuracy_data %>%
      dplyr::group_by(tool) %>%
      dplyr::arrange(desc(f1)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(f1))

    display_df <- data.frame(
      Tool = summary_df$tool,
      Window = summary_df$window,
      F1 = sprintf("%.3f", summary_df$f1),
      Precision = sprintf("%.3f", summary_df$precision),
      Recall = sprintf("%.3f", summary_df$recall),
      check.names = FALSE
    )

    if ("auprc" %in% names(summary_df)) {
      display_df$AUPRC <- sprintf("%.3f", summary_df$auprc)
    }
    if ("auroc" %in% names(summary_df)) {
      display_df$AUROC <- sprintf("%.3f", summary_df$auroc)
    }

    display_df
  }
}

# =============================================================================
# Figure 1: Summary Table as Plot
# =============================================================================

create_summary_table_plot <- function(optimal_data, accuracy_data) {
  summary_df <- create_optimal_summary_table(optimal_data, accuracy_data)

  # Convert to table grob
  table_grob <- gridExtra::tableGrob(
    summary_df,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = 8,
      padding = grid::unit(c(3, 3), "mm"),
      core = list(
        bg_params = list(fill = "white"),
        fg_params = list(fontsize = 8)
      ),
      colhead = list(
        bg_params = list(fill = "#2E86AB"),
        fg_params = list(fontface = "bold", col = "white", fontsize = 9)
      ),
      rowhead = list(
        fg_params = list(fontsize = 8)
      )
    )
  )

  # Add title
  title <- grid::textGrob(
    "Optimal Configuration Summary",
    x = grid::unit(0.5, "npc"),
    hjust = 0.5,
    gp = grid::gpar(fontsize = 10, fontface = "bold")
  )

  combined <- gridExtra::arrangeGrob(
    table_grob,
    top = title
  )

  combined
}

# =============================================================================
# Figure 2: Score Column Usage Summary
# =============================================================================

create_score_column_summary <- function(optimal_data) {
  if (is.null(optimal_data) || !"score_column" %in% names(optimal_data)) {
    message("Score column data not available")
    return(NULL)
  }

  # Count usage of each score column
  score_counts <- optimal_data %>%
    dplyr::mutate(
      score_label = ifelse(
        is.na(score_column) | score_column == "",
        "default",
        score_column
      )
    ) %>%
    dplyr::group_by(score_label) %>%
    dplyr::summarise(
      count = dplyr::n(),
      mean_f1 = mean(f1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(count))

  score_counts$score_label <- factor(score_counts$score_label,
                                      levels = score_counts$score_label)

  p <- ggplot(score_counts, aes(x = score_label, y = count, fill = mean_f1)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = count), vjust = -0.5, size = 2.5) +
    scale_fill_gradient(
      low = "#F0E442", high = "#009E73",
      name = "Mean F1"
    ) +
    labs(
      title = "Optimal Score Column Selection",
      subtitle = "How many tools perform best with each score column",
      x = "Score Column",
      y = "Number of Tools"
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(score_counts$count) * 1.2)) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

# =============================================================================
# Figure 3: Tool Performance Ranking
# =============================================================================

create_ranking_plot <- function(optimal_data, accuracy_data) {
  if (!is.null(optimal_data) && nrow(optimal_data) > 0) {
    ranking_df <- optimal_data %>%
      dplyr::select(tool, f1, precision, recall) %>%
      dplyr::arrange(desc(f1))
  } else {
    ranking_df <- accuracy_data %>%
      dplyr::group_by(tool) %>%
      dplyr::summarise(
        f1 = mean(f1, na.rm = TRUE),
        precision = mean(precision, na.rm = TRUE),
        recall = mean(recall, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(desc(f1))
  }

  ranking_df$tool <- factor(ranking_df$tool, levels = ranking_df$tool)
  ranking_df$rank <- seq_len(nrow(ranking_df))

  # Create ranking visualization
  df_long <- ranking_df %>%
    tidyr::pivot_longer(
      cols = c(f1, precision, recall),
      names_to = "metric",
      values_to = "value"
    )

  df_long$metric <- factor(df_long$metric,
                           levels = c("f1", "precision", "recall"),
                           labels = c("F1", "Precision", "Recall"))

  p <- ggplot(df_long, aes(x = tool, y = value, fill = metric)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(
      values = c("F1" = "#009E73", "Precision" = "#0072B2", "Recall" = "#D55E00")
    ) +
    geom_text(aes(label = sprintf("%.2f", value)),
              position = position_stack(vjust = 0.5),
              size = 2, color = "white") +
    labs(
      title = "Tool Performance Ranking",
      subtitle = "Tools ordered by F1 score (optimal configuration)",
      x = "Tool",
      y = "Score",
      fill = "Metric"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  p
}

# =============================================================================
# Figure 4: Metric Comparison Table
# =============================================================================

create_metric_comparison_table <- function(accuracy_data) {
  # Aggregate by tool
  summary_df <- accuracy_data %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(
      f1 = mean(f1, na.rm = TRUE),
      precision = mean(precision, na.rm = TRUE),
      recall = mean(recall, na.rm = TRUE),
      tp = sum(tp, na.rm = TRUE),
      fp = sum(fp, na.rm = TRUE),
      fn = sum(fn, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(f1)) %>%
    dplyr::mutate(
      f1 = sprintf("%.3f", f1),
      precision = sprintf("%.3f", precision),
      recall = sprintf("%.3f", recall)
    )

  display_df <- data.frame(
    Tool = summary_df$tool,
    F1 = summary_df$f1,
    Precision = summary_df$precision,
    Recall = summary_df$recall,
    TP = format(summary_df$tp, big.mark = ","),
    FP = format(summary_df$fp, big.mark = ","),
    FN = format(summary_df$fn, big.mark = ","),
    check.names = FALSE
  )

  # Create table grob
  table_grob <- gridExtra::tableGrob(
    display_df,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = 8,
      padding = grid::unit(c(2, 2), "mm"),
      core = list(
        bg_params = list(fill = "white"),
        fg_params = list(fontsize = 8)
      ),
      colhead = list(
        bg_params = list(fill = "#2E86AB"),
        fg_params = list(fontface = "bold", col = "white")
      )
    )
  )

  table_grob
}

# =============================================================================
# Figure 5: Best Tool per Metric Highlight
# =============================================================================

create_best_metric_highlight <- function(optimal_data, accuracy_data) {
  if (!is.null(optimal_data) && nrow(optimal_data) > 0) {
    base_df <- optimal_data
  } else {
    base_df <- accuracy_data %>%
      dplyr::group_by(tool) %>%
      dplyr::summarise(
        f1 = mean(f1, na.rm = TRUE),
        precision = mean(precision, na.rm = TRUE),
        recall = mean(recall, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # Find best for each metric
  metrics <- c("f1", "precision", "recall")
  best <- list()
  for (m in metrics) {
    best[[m]] <- base_df %>%
      dplyr::slice(which.max(get(m))) %>%
      dplyr::pull(tool)
  }

  # Create highlight data frame
  highlight_df <- data.frame(
    Metric = c("F1 Score", "Precision", "Recall"),
    Best_Tool = best,
    Score = c(
      max(base_df$f1, na.rm = TRUE),
      max(base_df$precision, na.rm = TRUE),
      max(base_df$recall, na.rm = TRUE)
    )
  )

  highlight_df$Score <- sprintf("%.3f", highlight_df$Score)

  # Create visual
  p <- ggplot(highlight_df, aes(x = Metric, y = as.numeric(Score), fill = Best_Tool)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = paste0(Best_Tool, "\n", Score)),
              vjust = 0.5, size = 2.5, color = "white") +
    scale_fill_manual(values = assign_tool_colors(unique(highlight_df$Best_Tool))) +
    labs(
      title = "Best Tool by Metric",
      subtitle = "Tool with highest score for each metric",
      x = "",
      y = "Score"
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1)) +
    theme_nature() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 9)
    )

  p
}

# =============================================================================
# Generate and Save Figures
# =============================================================================

message("Generating optimal metrics summary figures...")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Figure 1: Summary table
p1 <- create_summary_table_plot(optimal_df, accuracy_df)
if (!is.null(p1)) {
  save_figure(p1, file.path(output_dir, "summary_01_optimal_table.pdf"),
              width = 140, height = 100)
}

# Figure 2: Score column summary
p2 <- create_score_column_summary(optimal_df)
if (!is.null(p2)) {
  save_figure(p2, file.path(output_dir, "summary_02_score_columns.pdf"),
              width = 100, height = 70)
  save_figure(p2, file.path(output_dir, "summary_02_score_columns.png"),
              width = 100, height = 70, dpi = 300)
}

# Figure 3: Ranking plot
p3 <- create_ranking_plot(optimal_df, accuracy_df)
save_figure(p3, file.path(output_dir, "summary_03_ranking.pdf"),
            width = 120, height = 85)
save_figure(p3, file.path(output_dir, "summary_03_ranking.png"),
            width = 120, height = 85, dpi = 300)

# Figure 4: Comparison table
p4 <- create_metric_comparison_table(accuracy_df)
if (!is.null(p4)) {
  save_figure(p4, file.path(output_dir, "summary_04_comparison_table.pdf"),
              width = 140, height = 100)
}

# Figure 5: Best highlight
p5 <- create_best_metric_highlight(optimal_df, accuracy_df)
save_figure(p5, file.path(output_dir, "summary_05_best_by_metric.pdf"),
            width = 85, height = 85)
save_figure(p5, file.path(output_dir, "summary_05_best_by_metric.png"),
            width = 85, height = 85, dpi = 300)

message("Done! Optimal metrics summary figures saved to: ", output_dir)
