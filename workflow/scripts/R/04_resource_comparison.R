#!/usr/bin/env Rscript
# =============================================================================
# Resource Comparison Figure
# =============================================================================
# Generates computational resource comparison visualizations:
# - Runtime comparison (bar chart, log scale if needed)
# - Memory usage comparison
# - I/O comparison
# - Combined efficiency plot (F1 vs runtime)
#
# Usage:
#   Rscript 04_resource_comparison.R --resources resource_summary.tsv --accuracy accuracy_summary.tsv --output figures/
# =============================================================================

options(echo = FALSE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
resource_file <- "resource_summary.tsv"
accuracy_file <- "accuracy_summary.tsv"
output_dir <- "figures"

for (i in seq_along(args)) {
  if (args[i] == "--resources" && i + 1 <= length(args)) {
    resource_file <- args[i + 1]
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

message("Loading resource data from: ", resource_file)
resources <- load_resource_summary(resource_file)

message("Loading accuracy data from: ", accuracy_file)
accuracy <- load_accuracy_summary(accuracy_file)

# Assign tool colors
tool_colors <- assign_tool_colors(unique(resources$tool))

# =============================================================================
# Data Preparation
# =============================================================================

# Aggregate resource data by tool
resource_summary <- resources %>%
  dplyr::group_by(tool) %>%
  dplyr::summarise(
    mean_time_s = mean(s, na.rm = TRUE),
    sd_time_s = sd(s, na.rm = TRUE),
    mean_memory_mb = mean(max_rss, na.rm = TRUE),
    sd_memory_mb = sd(max_rss, na.rm = TRUE),
    mean_cpu_time_s = mean(cpu_time, na.rm = TRUE),
    total_io_mb = sum(io_in, na.rm = TRUE) + sum(io_out, na.rm = TRUE),
    n_samples = dplyr::n(),
    .groups = "drop"
  )

# Get accuracy data for efficiency plot
accuracy_summary <- accuracy %>%
  dplyr::group_by(tool) %>%
  dplyr::summarise(mean_f1 = mean(f1, na.rm = TRUE), .groups = "drop")

# Combine
combined <- resource_summary %>%
  dplyr::left_join(accuracy_summary, by = "tool") %>%
  dplyr::mutate(
    memory_gb = mean_memory_mb / 1024,
    time_min = mean_time_s / 60
  )

# =============================================================================
# Figure 1: Runtime Comparison
# =============================================================================

create_runtime_plot <- function(data, colors) {
  plot_data <- data %>%
    dplyr::arrange(desc(mean_time_s))

  plot_data$tool <- factor(plot_data$tool, levels = plot_data$tool)
  plot_colors <- colors[match(levels(plot_data$tool), names(colors))]

  p <- ggplot(plot_data, aes(x = tool, y = mean_time_s, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = plot_colors, guide = "none") +
    geom_errorbar(aes(ymin = mean_time_s - sd_time_s,
                      ymax = mean_time_s + sd_time_s),
                  width = 0.2, linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.1f", mean_time_s)),
              vjust = -0.5, size = 2) +
    labs(
      title = "Wall-Clock Runtime by Tool",
      subtitle = "Error bars represent standard deviation across samples",
      x = "Tool",
      y = "Time (seconds)"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

create_runtime_log_plot <- function(data, colors) {
  plot_data <- data %>%
    dplyr::arrange(desc(mean_time_s))

  plot_data$tool <- factor(plot_data$tool, levels = plot_data$tool)
  plot_colors <- colors[match(levels(plot_data$tool), names(colors))]

  p <- ggplot(plot_data, aes(x = tool, y = mean_time_s, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = plot_colors, guide = "none") +
    geom_text(aes(label = sprintf("%.0f", mean_time_s)),
              vjust = -0.5, size = 2) +
    labs(
      title = "Wall-Clock Runtime by Tool (Log Scale)",
      subtitle = "Log scale used to visualize wide range of runtimes",
      x = "Tool",
      y = "Time (seconds, log scale)"
    ) +
    scale_y_log10() +
    annotation_logticks(side = "l") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

# =============================================================================
# Figure 2: Memory Usage Comparison
# =============================================================================

create_memory_plot <- function(data, colors) {
  plot_data <- data %>%
    dplyr::arrange(desc(mean_memory_mb))

  plot_data$tool <- factor(plot_data$tool, levels = plot_data$tool)
  plot_colors <- colors[match(levels(plot_data$tool), names(colors))]

  p <- ggplot(plot_data, aes(x = tool, y = memory_gb, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = plot_colors, guide = "none") +
    geom_errorbar(aes(ymin = memory_gb - (sd_memory_mb / 1024),
                      ymax = memory_gb + (sd_memory_mb / 1024)),
                  width = 0.2, linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.1f", memory_gb)),
              vjust = -0.5, size = 2) +
    labs(
      title = "Peak Memory Usage by Tool",
      subtitle = "Error bars represent standard deviation across samples",
      x = "Tool",
      y = "Memory (GB)"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

# =============================================================================
# Figure 3: I/O Comparison
# =============================================================================

create_io_plot <- function(data, colors) {
  if (!"total_io_mb" %in% names(data) || all(is.na(data$total_io_mb))) {
    message("I/O data not available")
    return(NULL)
  }

  plot_data <- data %>%
    dplyr::arrange(desc(total_io_mb))

  plot_data$tool <- factor(plot_data$tool, levels = plot_data$tool)
  plot_colors <- colors[match(levels(plot_data$tool), names(colors))]

  p <- ggplot(plot_data, aes(x = tool, y = total_io_mb / 1024, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = plot_colors, guide = "none") +
    geom_text(aes(label = sprintf("%.1f", total_io_mb / 1024)),
              vjust = -0.5, size = 2) +
    labs(
      title = "Total I/O by Tool",
      subtitle = "Combined read and write operations",
      x = "Tool",
      y = "I/O (GB)"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

# =============================================================================
# Figure 4: Combined Efficiency (F1 vs Resources)
# =============================================================================

create_efficiency_scatter <- function(data, colors) {
  # Efficiency: F1 per unit of resource
  df_eff <- data %>%
    tidyr::pivot_longer(
      cols = c(time_min, memory_gb),
      names_to = "resource_type",
      values_to = "resource_value"
    ) %>%
    dplyr::filter(!is.na(resource_value), resource_value > 0)

  df_eff$resource_label <- ifelse(
    df_eff$resource_type == "time_min",
    "Time (min)",
    "Memory (GB)"
  )

  p <- ggplot(df_eff,
              aes(x = resource_value, y = mean_f1, color = tool, shape = tool)) +
    geom_point(size = 2.5) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = 1:length(unique(df_eff$tool))) +
    facet_wrap(~ resource_label, scales = "free_x") +
    labs(
      title = "Efficiency Analysis: F1 Score vs Resource Usage",
      subtitle = "Top-left corner = best performance with lowest resource use",
      x = "Resource Usage",
      y = "Mean F1 Score",
      color = "Tool",
      shape = "Tool"
    ) +
    theme_nature() +
    theme(
      legend.position = "right",
      strip.text = element_text(size = 8)
    )

  p
}

# =============================================================================
# Figure 5: Resource Profile Radar
# =============================================================================

create_resource_radar <- function(data) {
  # Normalize resources (0-1 scale, inverted so lower is better)
  metrics <- c("mean_time_s", "mean_memory_mb", "total_io_mb")
  available <- intersect(metrics, names(data))

  if (length(available) < 2) {
    message("Not enough resource metrics for radar plot")
    return(NULL)
  }

  df_norm <- data[, c("tool", available)]
  for (col in available) {
    max_val <- max(df_norm[[col]], na.rm = TRUE)
    if (max_val > 0) {
      # Invert: 1 - (value/max) so lower resource use = higher score
      df_norm[[col]] <- 1 - (df_norm[[col]] / max_val)
    }
  }

  # Add F1 (higher is better, already in correct direction)
  if ("mean_f1" %in% names(data)) {
    df_norm$mean_f1 <- data$mean_f1
  }

  df_long <- df_norm %>%
    tidyr::pivot_longer(
      cols = -tool,
      names_to = "metric",
      values_to = "value"
    )

  # Rename metrics for display
  df_long$metric <- factor(df_long$metric,
                          levels = c("mean_f1", "mean_time_s", "mean_memory_mb", "total_io_mb"),
                          labels = c("F1 Score", "Speed", "Memory Eff.", "I/O Eff."))

  p <- ggplot(df_long, aes(x = metric, y = value, fill = metric)) +
    geom_bar(stat = "identity", width = 0.7) +
    facet_wrap(~ tool) +
    scale_fill_manual(
      values = c("F1 Score" = "#009E73",
                 "Speed" = "#0072B2",
                 "Memory Eff." = "#56B4E9",
                 "I/O Eff." = "#E69F00"),
      guide = "none"
    ) +
    labs(
      title = "Resource Efficiency Profile",
      subtitle = "Higher scores indicate better efficiency",
      x = "",
      y = "Normalized Score (0-1)"
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_nature() +
    theme(
      strip.text = element_text(size = 7),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
    )

  p
}

# =============================================================================
# Figure 6: Comprehensive Resource Summary
# =============================================================================

create_resource_summary_multi <- function(data, colors) {
  plot_data <- data %>%
    dplyr::mutate(
      cpu_hours = mean_cpu_time_s / 3600
    )

  plot_data$tool <- factor(plot_data$tool,
                            levels = plot_data$tool[order(plot_data$mean_time_s,
                                                         decreasing = TRUE)])

  # Create a multi-panel plot
  p1 <- ggplot(plot_data, aes(x = tool, y = time_min, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = colors, guide = "none") +
    labs(title = "Runtime", x = "", y = "Time (min)") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.title.y = element_text(size = 7))

  p2 <- ggplot(plot_data, aes(x = tool, y = memory_gb, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = colors, guide = "none") +
    labs(title = "Memory", x = "", y = "Memory (GB)") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.title.y = element_text(size = 7))

  p3 <- ggplot(plot_data, aes(x = tool, y = cpu_hours, fill = tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = colors, guide = "none") +
    labs(title = "CPU Time", x = "", y = "CPU Hours") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.title.y = element_text(size = 7))

  # Combine
  combined <- gridExtra::grid.arrange(p1, p2, p3, nrow = 1)

  combined
}

# =============================================================================
# Generate and Save Figures
# =============================================================================

message("Generating resource comparison figures...")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Figure 1: Runtime
p1 <- create_runtime_plot(combined, tool_colors)
save_figure(p1, file.path(output_dir, "resource_01_runtime.pdf"),
            width = 100, height = 70)
save_figure(p1, file.path(output_dir, "resource_01_runtime.png"),
            width = 100, height = 70, dpi = 300)

# Figure 1b: Runtime (log scale)
p1b <- create_runtime_log_plot(combined, tool_colors)
save_figure(p1b, file.path(output_dir, "resource_01b_runtime_log.pdf"),
            width = 100, height = 70)

# Figure 2: Memory
p2 <- create_memory_plot(combined, tool_colors)
save_figure(p2, file.path(output_dir, "resource_02_memory.pdf"),
            width = 100, height = 70)
save_figure(p2, file.path(output_dir, "resource_02_memory.png"),
            width = 100, height = 70, dpi = 300)

# Figure 3: I/O
p3 <- create_io_plot(combined, tool_colors)
if (!is.null(p3)) {
  save_figure(p3, file.path(output_dir, "resource_03_io.pdf"),
              width = 100, height = 70)
  save_figure(p3, file.path(output_dir, "resource_03_io.png"),
              width = 100, height = 70, dpi = 300)
}

# Figure 4: Efficiency scatter
p4 <- create_efficiency_scatter(combined, tool_colors)
save_figure(p4, file.path(output_dir, "resource_04_efficiency.pdf"),
            width = 140, height = 85)
save_figure(p4, file.path(output_dir, "resource_04_efficiency.png"),
            width = 140, height = 85, dpi = 300)

# Figure 5: Resource profile
p5 <- create_resource_radar(combined)
if (!is.null(p5)) {
  save_figure(p5, file.path(output_dir, "resource_05_profile.pdf"),
              width = 140, height = 100)
  save_figure(p5, file.path(output_dir, "resource_05_profile.png"),
              width = 140, height = 100, dpi = 300)
}

# Figure 6: Multi-panel summary
p6 <- create_resource_summary_multi(combined, tool_colors)
if (!is.null(p6)) {
  save_figure(p6, file.path(output_dir, "resource_06_summary.pdf"),
              width = 174, height = 60)
}

message("Done! Resource comparison figures saved to: ", output_dir)
