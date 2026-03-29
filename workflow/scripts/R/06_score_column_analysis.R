#!/usr/bin/env Rscript
# =============================================================================
# Score Column Analysis Figures
# =============================================================================
# Generates per-tool score column comparison visualizations:
# - Fig A: Faceted grouped bar chart (AUROC/AUPRC/F1 per score column per tool)
# - Fig B: Best score selection dot plot (summary across all tools)
# - Fig C: Per-tool individual panels (detailed lollipop chart)
#
# Input: score_column_summary.tsv from benchmark_score_optimization.py
# =============================================================================

# =============================================================================
# Figure: Per-Tool Score Column Comparison (Faceted)
# =============================================================================
# Shows AUROC, AUPRC, and F1 for every candidate score column within each tool.
# Best score column is highlighted with a bold border and star annotation.

create_fig_score_columns_faceted <- function(summary_df, tool_colors) {
  if (is.null(summary_df) || nrow(summary_df) == 0) return(NULL)

  # Aggregate across modification types if present
  df_agg <- summary_df %>%
    dplyr::group_by(tool, score_column, score_type, is_best) %>%
    dplyr::summarise(
      auroc = mean(auroc, na.rm = TRUE),
      auprc = mean(auprc, na.rm = TRUE),
      f1 = mean(optimal_f1, na.rm = TRUE),
      .groups = "drop"
    )

  # Pivot to long format for grouped bars
  df_long <- df_agg %>%
    tidyr::pivot_longer(
      cols = c(auroc, auprc, f1),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = factor(metric,
                       levels = c("auroc", "auprc", "f1"),
                       labels = c("AUROC", "AUPRC", "F1"))
    )

  # Truncate long score column names for display
  df_long <- df_long %>%
    dplyr::mutate(
      score_label = ifelse(nchar(score_column) > 15,
                            paste0(substr(score_column, 1, 13), ".."),
                            score_column)
    )

  # Order score columns within each tool by AUROC (descending)
  score_order <- df_agg %>%
    dplyr::arrange(tool, desc(auroc)) %>%
    dplyr::group_by(tool) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(tool, score_column, rank)

  df_long <- df_long %>%
    dplyr::left_join(score_order, by = c("tool", "score_column")) %>%
    dplyr::mutate(
      score_label = forcats::fct_reorder2(score_label, tool, rank, .fun = function(t, r) -min(r))
    )

  # Reorder score_label within each tool panel by rank
  # Since facets share factor levels, we create a combined label
  df_long <- df_long %>%
    dplyr::arrange(tool, rank) %>%
    dplyr::mutate(
      panel_label = paste0(tool, "___", score_label),
      panel_label = forcats::fct_inorder(panel_label)
    )

  # Compute per-tool label ordering
  label_orders <- df_long %>%
    dplyr::distinct(tool, score_column, score_label, rank) %>%
    dplyr::arrange(tool, rank)

  # For facet_wrap, reorder score_label per tool
  # Use interaction trick: within each facet, order by rank
  df_long <- df_long %>%
    dplyr::group_by(tool) %>%
    dplyr::mutate(
      score_label = forcats::fct_reorder(score_label, rank, .fun = min)
    ) %>%
    dplyr::ungroup()

  # Color palette for metrics
  metric_colors <- c("AUROC" = "#0072B2", "AUPRC" = "#D55E00", "F1" = "#009E73")

  # Determine number of tools for layout
  n_tools <- length(unique(df_long$tool))
  ncol_facet <- min(3, n_tools)

  p <- ggplot(df_long, aes(x = score_label, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = NA) +
    # Highlight best score column with bold outline
    geom_bar(
      data = df_long %>% dplyr::filter(is_best == TRUE | is_best == "True"),
      stat = "identity", position = position_dodge(width = 0.8),
      width = 0.7, color = "black", linewidth = 0.6, fill = NA
    ) +
    # Add star for best
    geom_text(
      data = df_long %>%
        dplyr::filter((is_best == TRUE | is_best == "True") & metric == "AUROC") %>%
        dplyr::group_by(tool, score_column) %>%
        dplyr::slice(1),
      aes(x = score_label, y = 1.02, label = "\u2605"),
      inherit.aes = FALSE, size = 3, color = "#E69F00"
    ) +
    scale_fill_manual(values = metric_colors, name = "Metric") +
    facet_wrap(~tool, scales = "free_x", ncol = ncol_facet) +
    labs(
      title = "Score Column Performance by Tool",
      subtitle = "AUROC, AUPRC, and F1 for each candidate score column (\u2605 = best)",
      x = "Score Column",
      y = "Score"
    ) +
    scale_y_continuous(limits = c(0, 1.08), expand = c(0, 0)) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1, size = 6),
      legend.position = "top",
      strip.text = element_text(size = 7, face = "bold"),
      panel.spacing = unit(3, "mm")
    )

  p
}


# =============================================================================
# Figure: Best Score Selection Summary (Dot Plot)
# =============================================================================
# For each tool, shows all score columns as dots (AUROC on x-axis),
# with the best score column highlighted as a larger filled dot.

create_fig_best_score_dotplot <- function(summary_df, tool_colors) {
  if (is.null(summary_df) || nrow(summary_df) == 0) return(NULL)

  # Aggregate across modification types
  df_agg <- summary_df %>%
    dplyr::group_by(tool, score_column, score_type, is_best) %>%
    dplyr::summarise(
      auroc = mean(auroc, na.rm = TRUE),
      auprc = mean(auprc, na.rm = TRUE),
      f1 = mean(optimal_f1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(auroc))

  if (nrow(df_agg) == 0) return(NULL)

  # Order tools by best AUROC
  tool_order <- df_agg %>%
    dplyr::filter(is_best == TRUE | is_best == "True") %>%
    dplyr::arrange(desc(auroc)) %>%
    dplyr::pull(tool) %>%
    unique()
  # Add any tools without a best flag
  remaining <- setdiff(unique(df_agg$tool), tool_order)
  tool_order <- c(tool_order, remaining)

  df_agg$tool <- factor(df_agg$tool, levels = rev(tool_order))

  # Truncate score column names
  df_agg <- df_agg %>%
    dplyr::mutate(
      label = ifelse(nchar(score_column) > 18,
                      paste0(substr(score_column, 1, 16), ".."),
                      score_column),
      is_best_flag = (is_best == TRUE | is_best == "True")
    )

  p <- ggplot(df_agg, aes(x = auroc, y = tool)) +
    # Non-best points: small, gray
    geom_point(
      data = df_agg %>% dplyr::filter(!is_best_flag),
      aes(shape = score_type),
      size = 2, color = "gray60", alpha = 0.6
    ) +
    # Best points: large, colored
    geom_point(
      data = df_agg %>% dplyr::filter(is_best_flag),
      aes(color = tool, shape = score_type),
      size = 4, stroke = 0.8
    ) +
    # Label best score column
    ggrepel::geom_text_repel(
      data = df_agg %>% dplyr::filter(is_best_flag),
      aes(label = label),
      size = 2.5, nudge_x = 0.02, max.overlaps = 20,
      segment.size = 0.2, segment.color = "gray50"
    ) +
    # Label non-best with smaller text
    ggrepel::geom_text_repel(
      data = df_agg %>% dplyr::filter(!is_best_flag),
      aes(label = label),
      size = 1.8, color = "gray50", nudge_x = 0.01,
      max.overlaps = 15, segment.size = 0.15, segment.color = "gray70"
    ) +
    scale_color_manual(values = tool_colors, guide = "none") +
    scale_shape_manual(
      values = c("p-value" = 17, "probability" = 16),
      name = "Score Type"
    ) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray80", linewidth = 0.25) +
    labs(
      title = "Best Score Column Selection",
      subtitle = "AUROC for all candidate score columns per tool (large dot = selected best)",
      x = "AUROC",
      y = "Tool"
    ) +
    scale_x_continuous(limits = c(max(0.4, min(df_agg$auroc, na.rm = TRUE) - 0.05), 1.0)) +
    theme_nature() +
    theme(
      legend.position = "top",
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.15)
    )

  p
}


# =============================================================================
# Figure: Per-Tool Score Column Heatmap
# =============================================================================
# Heatmap with tools on y-axis, score columns on x-axis, fill = AUROC.
# Best score column annotated.

create_fig_score_heatmap <- function(summary_df, tool_colors) {
  if (is.null(summary_df) || nrow(summary_df) == 0) return(NULL)

  # Aggregate
  df_agg <- summary_df %>%
    dplyr::group_by(tool, score_column, is_best) %>%
    dplyr::summarise(
      auroc = mean(auroc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(auroc))

  if (nrow(df_agg) == 0) return(NULL)

  # Order tools by best AUROC
  tool_order <- df_agg %>%
    dplyr::filter(is_best == TRUE | is_best == "True") %>%
    dplyr::arrange(desc(auroc)) %>%
    dplyr::pull(tool) %>%
    unique()
  remaining <- setdiff(unique(df_agg$tool), tool_order)
  tool_order <- c(tool_order, remaining)
  df_agg$tool <- factor(df_agg$tool, levels = rev(tool_order))

  # Truncate column names
  df_agg <- df_agg %>%
    dplyr::mutate(
      score_label = ifelse(nchar(score_column) > 20,
                            paste0(substr(score_column, 1, 18), ".."),
                            score_column),
      is_best_flag = (is_best == TRUE | is_best == "True")
    )

  p <- ggplot(df_agg, aes(x = score_label, y = tool, fill = auroc)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", auroc)),
              size = 2, color = "black") +
    # Bold border for best
    geom_tile(
      data = df_agg %>% dplyr::filter(is_best_flag),
      color = "black", linewidth = 0.8, fill = NA
    ) +
    scale_fill_viridis_c(
      option = "viridis", limits = c(0.4, 1.0), na.value = "gray90",
      name = "AUROC"
    ) +
    labs(
      title = "Score Column AUROC Heatmap",
      subtitle = "AUROC per tool \u00d7 score column (bold border = best)",
      x = "Score Column",
      y = "Tool"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1, size = 6),
      legend.position = "right",
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank()
    )

  p
}


# =============================================================================
# Figure: Per-Tool Lollipop Charts (Individual Panels)
# =============================================================================
# For each tool, a lollipop chart showing score columns ranked by AUROC,
# with F1 annotated. This gives the most detailed per-tool view.

create_fig_score_lollipop_per_tool <- function(summary_df, tool_colors) {
  if (is.null(summary_df) || nrow(summary_df) == 0) return(NULL)

  # Aggregate
  df_agg <- summary_df %>%
    dplyr::group_by(tool, score_column, score_type, is_best) %>%
    dplyr::summarise(
      auroc = mean(auroc, na.rm = TRUE),
      auprc = mean(auprc, na.rm = TRUE),
      f1 = mean(optimal_f1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(auroc))

  if (nrow(df_agg) == 0) return(NULL)

  tools <- sort(unique(df_agg$tool))
  plot_list <- list()

  for (t in tools) {
    df_tool <- df_agg %>%
      dplyr::filter(tool == t) %>%
      dplyr::arrange(desc(auroc)) %>%
      dplyr::mutate(
        is_best_flag = (is_best == TRUE | is_best == "True"),
        score_label = ifelse(nchar(score_column) > 22,
                              paste0(substr(score_column, 1, 20), ".."),
                              score_column),
        score_label = forcats::fct_inorder(score_label)
      )

    # Skip tools with only 1 score column
    if (nrow(df_tool) < 2) next

    tool_color <- if (t %in% names(tool_colors)) tool_colors[t] else "#0072B2"

    p <- ggplot(df_tool, aes(x = auroc, y = forcats::fct_rev(score_label))) +
      # Lollipop stem
      geom_segment(
        aes(x = 0.5, xend = auroc, yend = forcats::fct_rev(score_label)),
        color = "gray70", linewidth = 0.4
      ) +
      # AUROC points
      geom_point(
        aes(size = f1, fill = is_best_flag),
        shape = 21, stroke = 0.5, color = "black"
      ) +
      # F1 annotation
      geom_text(
        aes(label = sprintf("F1=%.2f", f1)),
        hjust = -0.3, size = 2, color = "gray30"
      ) +
      scale_fill_manual(
        values = c("FALSE" = "gray80", "TRUE" = tool_color),
        guide = "none"
      ) +
      scale_size_continuous(range = c(2, 5), name = "F1", guide = "none") +
      labs(
        title = t,
        x = "AUROC",
        y = NULL
      ) +
      scale_x_continuous(
        limits = c(0.45, min(1.15, max(df_tool$auroc, na.rm = TRUE) + 0.15)),
        expand = c(0, 0)
      ) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray80", linewidth = 0.2) +
      theme_nature() +
      theme(
        plot.title = element_text(size = 8, face = "bold", color = tool_color),
        axis.text.y = element_text(size = 6),
        plot.margin = margin(2, 8, 2, 2, "mm")
      )

    plot_list[[t]] <- p
  }

  if (length(plot_list) == 0) return(NULL)

  # Arrange in a grid
  n <- length(plot_list)
  ncol_grid <- min(3, n)
  nrow_grid <- ceiling(n / ncol_grid)

  combined <- gridExtra::arrangeGrob(
    grobs = plot_list,
    ncol = ncol_grid,
    top = grid::textGrob(
      "Score Column Ranking by Tool (AUROC, dot size = F1)",
      gp = grid::gpar(fontsize = 9, fontface = "bold")
    )
  )

  combined
}
