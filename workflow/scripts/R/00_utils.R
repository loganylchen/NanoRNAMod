# =============================================================================
# Nature-Style R Visualization Utilities
# =============================================================================
# Shared utility functions for Nature-quality figures
#
# Usage:
#   source("00_utils.R")
#   plot <- ggplot(...) + theme_nature()
#   save_figure(plot, "figure1.pdf", width = 85)
# =============================================================================

required_packages <- c("ggplot2", "grid", "gridExtra", "scales",
                       "RColorBrewer", "viridis", "readr", "dplyr",
                       "tidyr", "purrr")

install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("Installing package:", pkg))
      install.packages(pkg, repos = "https://cloud.r-project.org/")
      library(pkg, character.only = TRUE)
    }
  }
}

install_if_missing(required_packages)

# =============================================================================
# Color Palettes (Colorblind-Friendly)
# =============================================================================

#' Okabe-Ito Color Palette
#' Colorblind-friendly palette from Okabe & Ito (2008)
#' @param n Number of colors to return
#' @return Character vector of hex colors
okabe_ito_palette <- function(n = NULL) {
  colors <- c(
    "black" = "#000000",
    "orange" = "#E69F00",
    "sky_blue" = "#56B4E9",
    "bluish_green" = "#009E73",
    "yellow" = "#F0E442",
    "blue" = "#0072B2",
    "vermilion" = "#D55E00",
    "reddish_purple" = "#CC79A7",
    "gray" = "#999999"
  )
  if (is.null(n)) return(colors)
  if (n > length(colors)) {
    warning(paste("Requested", n, "colors but only", length(colors), "available. Recycling."))
  }
  rep(colors, length.out = n)
}

#' Nature Publication Colors
#' Common colors used in Nature figures
#' @return List of color vectors
nature_colors <- list(
  primary = okabe_ito_palette()[c("blue", "vermilion", "bluish_green", "orange")],
  sequential = colorRampPalette(c("#F0E442", "#F39C12", "#E74C3C", "#8E44AD")),
  diverging = colorRampPalette(c("#2E86AB", "#FFFFFF", "#A23B72")),
  heatmap = colorRampPalette(c("#2E86AB", "#A8DADC", "#F1FAEE", "#F8B500", "#E63946"))
)

#' Tool-specific color assignment
#' @param tools Vector of tool names
#' @return Named character vector of colors
assign_tool_colors <- function(tools) {
  n_tools <- length(tools)
  colors <- okabe_ito_palette(n_tools)
  names(colors) <- tools
  # Sort tools alphabetically for consistent coloring
  colors[order(names(colors))]
}

# =============================================================================
# Nature Theme
# =============================================================================

#' Nature-style theme for ggplot2
#' Based on theme_classic() with Nature journal specifications
#' @param base_size Base font size (default 8 for Nature figures)
#' @param base_family Font family (Helvetica or Arial)
#' @return ggplot2 theme object
theme_nature <- function(base_size = 8, base_family = "Helvetica") {
  # Check if font is available, fallback to sans-serif
  # Handle case when systemfonts is not available
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    available_fonts <- systemfonts::system_fonts()$family
    if (!base_family %in% available_fonts) {
      if ("Arial" %in% available_fonts) {
        base_family <- "Arial"
      } else {
        base_family <- "sans"
      }
    }
  } else {
    # Fallback when systemfonts is not available
    base_family <- "sans"
  }

  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # Axis elements
      axis.line = element_line(color = "black", linewidth = 0.25, linetype = "solid"),
      axis.ticks = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length = unit(2.5, "pt"),
      axis.text = element_text(color = "black", size = base_size),

      # Panel and plot area
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25),

      # Legend
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.position = "right",

      # Grid (typically removed for Nature)
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # Strips
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = base_size, face = "bold"),

      # Plot margins (tight for Nature)
      plot.margin = margin(2, 2, 2, 2, "mm"),

      # Title
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0)
    )
}

#' Alternative theme with subtle grid (for some plot types)
theme_nature_grid <- function(base_size = 8, base_family = "Helvetica") {
  theme_nature(base_size, base_family) +
    theme(
      panel.grid.major = element_line(color = "gray80", linewidth = 0.1),
      panel.grid.minor = element_blank()
    )
}

# =============================================================================
# Save Functions
# =============================================================================

#' Save figure at publication resolution
#' @param plot ggplot object to save
#' @param filename Output filename (extension determines format)
#' @param width Width in mm (Nature single column: 85mm, double: 174mm)
#' @param height Height in mm
#' @param dpi Resolution (default 300 for publication)
#' @param ... Additional arguments to ggsave
save_figure <- function(plot, filename, width = 85, height = 85, dpi = 300, ...) {
  # Create output directory if needed
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)

  # Default to PDF if no extension
  if (!grepl("\\.(pdf|png|svg|eps|tiff)$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".pdf")
  }

  # Save with appropriate settings
  ggsave(
    filename = filename,
    plot = plot,
    width = width / 25.4,  # Convert mm to inches
    height = height / 25.4,
    dpi = dpi,
    units = "in",
    ...
  )

  message(paste("Saved:", filename))
  invisible(filename)
}

#' Save multi-panel figure
#' @param plots List of ggplot objects
#' @param filename Output filename
#' @param ncol Number of columns
#' @param nrow Number of rows
#' @param labels Panel labels (e.g., c("a", "b", "c"))
#' @param ... Arguments passed to save_figure
save_multipanel <- function(plots, filename, ncol = NULL, nrow = NULL,
                            labels = NULL, ...) {

  # Calculate grid dimensions
  n <- length(plots)
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- ceiling(sqrt(n))
    nrow <- ceiling(n / ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(n / nrow)
  } else if (is.null(nrow)) {
    nrow <- ceiling(n / ncol)
  }

  # Arrange plots
  combined <- gridExtra::grid.arrange(grobs = plots, ncol = ncol, nrow = nrow)

  # Add labels if provided
  if (!is.null(labels)) {
    combined <- gridExtra::arrangeGrob(
      combined,
      top = grid::textGrob(
        paste(labels, collapse = "  "),
        x = grid::unit(0, "npc"),
        hjust = 0,
        gp = grid::gpar(fontsize = 10, fontface = "bold")
      )
    )
  }

  # Calculate approximate dimensions
  width_mm <- 85 * ncol
  height_mm <- 85 * nrow

  save_figure(combined, filename, width = width_mm, height = height_mm, ...)
}

# =============================================================================
# Data Loading Functions
# =============================================================================

#' Load accuracy summary TSV
#' @param path Path to accuracy_summary.tsv or accuracy_summary_overall.tsv
#' @return Data frame with accuracy metrics
load_accuracy_summary <- function(path) {
  if (!file.exists(path)) {
    stop(paste("File not found:", path))
  }
  df <- readr::read_tsv(path, show_col_types = FALSE)

  # Ensure required columns exist
  required_cols <- c("tool", "f1", "precision", "recall")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    warning(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }

  df
}

#' Load optimal thresholds TSV
#' @param path Path to optimal_thresholds.tsv
#' @return Data frame with optimal thresholds
load_optimal_thresholds <- function(path) {
  if (!file.exists(path)) {
    stop(paste("File not found:", path))
  }
  readr::read_tsv(path, show_col_types = FALSE)
}

#' Load resource summary TSV
#' @param path Path to resource_summary.tsv
#' @return Data frame with resource metrics
load_resource_summary <- function(path) {
  if (!file.exists(path)) {
    stop(paste("File not found:", path))
  }
  df <- readr::read_tsv(path, show_col_types = FALSE)

  # Convert memory from MB to GB if needed
  if ("max_rss" %in% names(df)) {
    # Keep in MB for consistency
  }

  df
}

#' Load threshold evaluation TSV
#' @param path Path to threshold_evaluation.tsv
#' @return Data frame with threshold sweep results
load_threshold_evaluation <- function(path) {
  if (!file.exists(path)) {
    stop(paste("File not found:", path))
  }
  readr::read_tsv(path, show_col_types = FALSE)
}

#' Load all benchmark data from a directory
#' @param dir_path Directory containing benchmark TSV files
#' @return List with data frames: accuracy, overall, resources, thresholds, optimal
load_all_benchmark_data <- function(dir_path) {

  files <- list(
    accuracy = file.path(dir_path, "accuracy_summary.tsv"),
    overall = file.path(dir_path, "accuracy_summary_overall.tsv"),
    resources = file.path(dir_path, "resource_summary.tsv"),
    thresholds = file.path(dir_path, "threshold_evaluation.tsv"),
    optimal = file.path(dir_path, "optimal_thresholds.tsv")
  )

  result <- list()
  for (name in names(files)) {
    path <- files[[name]]
    if (file.exists(path)) {
      tryCatch({
        result[[name]] <- switch(name,
          accuracy = load_accuracy_summary(path),
          overall = load_accuracy_summary(path),
          resources = load_resource_summary(path),
          thresholds = load_threshold_evaluation(path),
          optimal = load_optimal_thresholds(path)
        )
        message(paste("Loaded:", name))
      }, error = function(e) {
        warning(paste("Could not load", name, ":", e$message))
      })
    } else {
      warning(paste("File not found:", path))
    }
  }

  result
}

# =============================================================================
# Data Transformation Helpers
# =============================================================================

#' Reshape metrics data for plotting
#' @param df Accuracy summary data frame
#' @return Reshaped data frame in long format
reshape_metrics_long <- function(df) {
  metrics <- c("precision", "recall", "f1", "auprc", "auroc", "mcc")
  available_metrics <- intersect(metrics, names(df))

  df %>%
    dplyr::select(tool, modification_type, window, any_of(available_metrics)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(available_metrics),
      names_to = "metric",
      values_to = "value"
    )
}

#' Calculate summary statistics by tool
#' @param df Accuracy summary data frame
#' @return Data frame with tool-level summaries
summarize_by_tool <- function(df) {
  df %>%
    dplyr::group_by(tool) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_f1 = mean(f1, na.rm = TRUE),
      sd_f1 = sd(f1, na.rm = TRUE),
      mean_precision = mean(precision, na.rm = TRUE),
      mean_recall = mean(recall, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(mean_f1))
}

#' Filter to optimal window per tool
#' @param df Accuracy summary data frame
#' @return Data frame filtered to best window for each tool
filter_optimal_window <- function(df) {
  df %>%
    dplyr::group_by(tool) %>%
    dplyr::arrange(desc(f1)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}

# =============================================================================
# Annotation Helpers
# =============================================================================

#' Add significance annotations to plot
#' @param plot ggplot object
#' @param comparisons List of comparisons to annotate
#' @param test_data Data frame with p-values
#' @param y_position Y position for annotations
#' @return ggplot with annotations
add_significance <- function(plot, comparisons, test_data, y_position = NULL) {
  # Implementation for adding significance stars
  # This is a placeholder for future implementation
  plot
}

#' Add optimal threshold annotation
#' @param plot ggplot object
#' @param threshold_value Optimal threshold value
#' @param label Annotation label
#' @return ggplot with annotation
add_threshold_annotation <- function(plot, threshold_value,
                                       label = "Optimal") {
  plot +
    geom_vline(xintercept = threshold_value,
               linetype = "dashed", color = "gray50",
               linewidth = 0.25) +
    annotate("text", x = threshold_value, y = Inf,
             label = label, vjust = 2, hjust = -0.1,
             size = 2.5, color = "gray30")
}

# =============================================================================
# Table Formatting
# =============================================================================

#' Format data frame as table grob for ggplot annotation
#' @param df Data frame to format
#' @param digits Numeric rounding
#' @return tableGrob object
format_table_grob <- function(df, digits = 3) {
  # Round numeric columns
  df_formatted <- df
  numeric_cols <- sapply(df, is.numeric)
  df_formatted[numeric_cols] <- lapply(df[numeric_cols], round, digits)

  # Create table grob
  gridExtra::tableGrob(
    df_formatted,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = 8,
      padding = grid::unit(c(2, 2), "mm"),
      core = list(bg_params = list(fill = "white")),
      colhead = list(bg_params = list(fill = "gray90"),
                     fg_params = list(fontface = "bold"))
    )
  )
}

# =============================================================================
# Export visible objects
# =============================================================================

# Export key functions for easy access
utils_export <- c(
  "theme_nature", "theme_nature_grid",
  "save_figure", "save_multipanel",
  "load_accuracy_summary", "load_optimal_thresholds",
  "load_resource_summary", "load_threshold_evaluation",
  "load_all_benchmark_data",
  "reshape_metrics_long", "summarize_by_tool", "filter_optimal_window",
  "assign_tool_colors", "okabe_ito_palette", "nature_colors",
  "add_threshold_annotation", "format_table_grob"
)

message("Nature visualization utilities loaded successfully.")
message(paste("Available functions:", paste(utils_export, collapse = ", ")))
