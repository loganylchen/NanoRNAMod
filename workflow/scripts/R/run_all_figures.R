#!/usr/bin/env Rscript
# =============================================================================
# Master Script: Generate All Benchmark Figures
# =============================================================================
# This script sources all individual figure scripts and generates
# a complete set of Nature-quality benchmark visualizations.
#
# Usage:
#   Rscript run_all_figures.R --data-dir /path/to/benchmark/data --output figures/
#   OR via Snakemake script directive (uses snakemake@input object)
#
# Expected files in data directory:
#   - accuracy_summary.tsv
#   - accuracy_summary_overall.tsv (optional)
#   - resource_summary.tsv
#   - threshold_evaluation.tsv
#   - optimal_thresholds.tsv
# =============================================================================

# Suppress startup messages
options(echo = FALSE, warn = 1)

# =============================================================================
# Configuration - Handle both Snakemake script and command-line invocation
# =============================================================================

# Get script directory - use alternative method that doesn't require this.file package
script_dir <- tryCatch({
  # Try to get from script path if available
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    path <- sub("^--file=", "", file_arg[1])
    dirname(normalizePath(path))
  } else {
    "."
  }
}, error = function(e) {
  # Fallback to current working directory
  "."
})

# Check if running as Snakemake script
if (exists("snakemake")) {
  # Snakemake script directive - use snakemake object
  input_files <- list(
    accuracy = snakemake@input[["accuracy"]],
    detailed = snakemake@input[["accuracy_overall"]],
    overall = snakemake@input[["accuracy_overall"]],
    optimal = snakemake@input[["optimal_scores"]],
    thresholds = snakemake@input[["thresholds"]],
    resources = snakemake@input[["resources"]]
  )
  output_dir <- snakemake@output[["dir"]]
  window_filter <- snakemake@params[["window"]]
  data_dir <- dirname(input_files$accuracy)
} else {
  # Command-line invocation
  data_dir <- "."
  output_dir <- "figures"
  window_filter <- NULL

  # Initialize input_files list
  input_files <- list(
    accuracy = NULL,
    detailed = NULL,
    overall = NULL,
    optimal = NULL,
    thresholds = NULL,
    resources = NULL
  )

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
    } else if (args[i] == "--accuracy" && i + 1 <= length(args)) {
      input_files$accuracy <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--accuracy-detailed" && i + 1 <= length(args)) {
      input_files$detailed <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--accuracy-overall" && i + 1 <= length(args)) {
      input_files$overall <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--optimal-scores" && i + 1 <= length(args)) {
      input_files$optimal <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--thresholds" && i + 1 <= length(args)) {
      input_files$thresholds <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--resources" && i + 1 <= length(args)) {
      input_files$resources <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  # Set default paths based on data_dir if not specified
  if (is.null(input_files$accuracy)) {
    input_files$accuracy <- file.path(data_dir, "accuracy_summary.tsv")
  }
  if (is.null(input_files$detailed)) {
    input_files$detailed <- file.path(data_dir, "accuracy_summary_overall.tsv")
  }
  if (is.null(input_files$overall)) {
    input_files$overall <- input_files$detailed
  }
  if (is.null(input_files$optimal)) {
    input_files$optimal <- file.path(data_dir, "optimal_score_per_tool.tsv")
  }
  if (is.null(input_files$thresholds)) {
    input_files$thresholds <- file.path(data_dir, "threshold_evaluation.tsv")
  }
  if (is.null(input_files$resources)) {
    input_files$resources <- file.path(data_dir, "resource_summary.tsv")
  }
}

# Define output subdirectories
output_dirs <- list(
  overview = file.path(output_dir, "01_performance_overview"),
  modification = file.path(output_dir, "02_per_modification"),
  threshold = file.path(output_dir, "03_threshold_optimization"),
  resource = file.path(output_dir, "04_resource_comparison"),
  summary = file.path(output_dir, "05_optimal_summary")
)

# =============================================================================
# Helper Functions
# =============================================================================

check_file_exists <- function(path, required = TRUE) {
  if (file.exists(path)) {
    return(TRUE)
  } else {
    if (required) {
      warning("Required file not found: ", path)
    }
    return(FALSE)
  }
}

run_script <- function(script_name, args, description) {
  message(paste("\n", paste(rep("=", 60), collapse = "")))
  message(paste("Running:", description))
  message(paste(rep("=", 60), collapse = ""))

  script_path <- file.path(script_dir, script_name)

  if (!file.exists(script_path)) {
    warning("Script not found: ", script_path)
    return(FALSE)
  }

  # Build command
  cmd <- paste("Rscript", script_path, args)

  # Run script
  result <- tryCatch({
    system(cmd, intern = FALSE)
    TRUE
  }, error = function(e) {
    warning("Error running script: ", e$message)
    FALSE
  })

  if (result) {
    message("Completed: ", description)
  }

  result
}

# =============================================================================
# Main Execution
# =============================================================================

message(paste(rep("=", 60), collapse = ""))
message("NanoRNAMod Benchmark Figure Generation")
message(paste(rep("=", 60), collapse = ""))

message("\nConfiguration:")
message("  Data directory: ", data_dir)
message("  Output directory: ", output_dir)
if (!is.null(window_filter)) {
  message("  Window filter: ", window_filter)
}

# Check input files
message("\nChecking input files...")
files_ok <- TRUE

if (!check_file_exists(input_files$accuracy, required = TRUE)) {
  files_ok <- FALSE
}

if (!file.exists(input_files$overall)) {
  message("  Note: accuracy_summary_overall.tsv not found, using accuracy_summary.tsv")
  input_files$overall <- input_files$accuracy
}

if (!check_file_exists(input_files$resources, required = FALSE)) {
  message("  Warning: resource_summary.tsv not found. Resource figures will be skipped.")
}

if (!check_file_exists(input_files$thresholds, required = FALSE)) {
  message("  Note: threshold_evaluation.tsv not found. Threshold figures will be limited.")
}

if (!check_file_exists(input_files$optimal, required = FALSE)) {
  message("  Note: optimal_thresholds.tsv not found. Using accuracy data for optimization.")
}

if (!files_ok) {
  stop("Required input files missing. Exiting.")
}

# Create output directories
message("\nCreating output directories...")
for (dir in output_dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  message("  Created: ", dir)
}

# =============================================================================
# Run Individual Scripts
# =============================================================================

results <- list()

# Script 1: Performance Overview
args_overview <- paste0(
  "--input ", input_files$accuracy,
  " --output ", output_dirs$overview
)
if (!is.null(window_filter)) {
  args_overview <- paste(args_overview, " --window ", window_filter)
}
results$overview <- run_script(
  "01_performance_overview.R",
  args_overview,
  "Performance Overview"
)

# Script 2: Per-Modification Analysis
args_mod <- paste0(
  "--input ", input_files$accuracy,
  " --output ", output_dirs$modification
)
results$modification <- run_script(
  "02_per_modification_analysis.R",
  args_mod,
  "Per-Modification Analysis"
)

# Script 3: Threshold Optimization (if data available)
if (file.exists(input_files$thresholds)) {
  args_thresh <- paste0(
    "--input ", input_files$thresholds,
    " --output ", output_dirs$threshold
  )
  results$threshold <- run_script(
    "03_threshold_optimization.R",
    args_thresh,
    "Threshold Optimization"
  )
} else {
  message("\nSkipping threshold optimization (no threshold_evaluation.tsv)")
  results$threshold <- NA
}

# Script 4: Resource Comparison (if data available)
if (file.exists(input_files$resources)) {
  args_res <- paste0(
    "--resources ", input_files$resources,
    " --accuracy ", input_files$overall,
    " --output ", output_dirs$resource
  )
  results$resource <- run_script(
    "04_resource_comparison.R",
    args_res,
    "Resource Comparison"
  )
} else {
  message("\nSkipping resource comparison (no resource_summary.tsv)")
  results$resource <- NA
}

# Script 5: Optimal Metrics Summary
args_summary <- paste0(
  "--optimal ", input_files$optimal,
  " --accuracy ", input_files$overall,
  " --output ", output_dirs$summary
)
# Handle case where optimal file doesn't exist
if (!file.exists(input_files$optimal)) {
  args_summary <- paste0(
    "--accuracy ", input_files$overall,
    " --output ", output_dirs$summary
  )
}
results$summary <- run_script(
  "05_optimal_metrics_summary.R",
  args_summary,
  "Optimal Metrics Summary"
)

# =============================================================================
# Summary Report
# =============================================================================

message("\n")
message(paste(rep("=", 60), collapse = ""))
message("Generation Summary")
message(paste(rep("=", 60), collapse = ""))

completed <- sum(!is.na(results) & results == TRUE)
total <- length(results)

message(sprintf("Scripts completed: %d/%d", completed, total))

if (results$overview) {
  message("  [OK] Performance Overview")
} else {
  message("  [SKIP] Performance Overview")
}

if (results$modification) {
  message("  [OK] Per-Modification Analysis")
} else {
  message("  [SKIP] Per-Modification Analysis")
}

if (!is.na(results$threshold)) {
  if (results$threshold) {
    message("  [OK] Threshold Optimization")
  } else {
    message("  [FAIL] Threshold Optimization")
  }
} else {
  message("  [SKIP] Threshold Optimization (no data)")
}

if (!is.na(results$resource)) {
  if (results$resource) {
    message("  [OK] Resource Comparison")
  } else {
    message("  [FAIL] Resource Comparison")
  }
} else {
  message("  [SKIP] Resource Comparison (no data)")
}

if (results$summary) {
  message("  [OK] Optimal Metrics Summary")
} else {
  message("  [SKIP] Optimal Metrics Summary")
}

message("\nOutput directory: ", output_dir)
message("\nFigure files:")
for (subdir in names(output_dirs)) {
  dir_path <- output_dirs[[subdir]]
  if (dir.exists(dir_path)) {
    files <- list.files(dir_path, pattern = "\\.(pdf|png)$")
    if (length(files) > 0) {
      message(sprintf("  %s: %d files", subdir, length(files)))
    }
  }
}

message("\n")
message(paste(rep("=", 60), collapse = ""))
message("Done!")
message(paste(rep("=", 60), collapse = ""))
