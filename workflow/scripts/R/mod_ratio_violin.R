#!/usr/bin/env Rscript
# Violin plot of modification ratio distributions
# Input: mod_ratio_data.tsv (from benchmark_mod_ratio.py)
# Output: mod_ratio_violin.pdf

library(ggplot2)
library(dplyr)

# Snakemake or CLI
if (exists("snakemake")) {
  input_data <- snakemake@input[["data"]]
  output_pdf <- snakemake@output[["pdf"]]
  log_file <- snakemake@log[[1]]
  log_con <- file(log_file, open = "wt")
  sink(log_con, type = "output")
  sink(log_con, type = "message")
} else {
  args <- commandArgs(trailingOnly = TRUE)
  input_data <- args[1]
  output_pdf <- args[2]
}

df <- read.delim(input_data, stringsAsFactors = FALSE)

if (nrow(df) == 0) {
  pdf(output_pdf, width = 8, height = 6)
  plot.new()
  text(0.5, 0.5, "No modification ratio data available", cex = 1.5)
  dev.off()
  quit(save = "no")
}

# Ensure mod_ratio is numeric and clamp to [0, 1]
df$mod_ratio <- as.numeric(df$mod_ratio)
df <- df[!is.na(df$mod_ratio), ]
df$mod_ratio <- pmin(pmax(df$mod_ratio, 0), 1)

# Count sites per group for x-axis labels
df <- df %>%
  group_by(tool, mode, truth) %>%
  mutate(
    n = n(),
    group_label = paste0(truth, "\n(n=", format(n, big.mark = ","), ")")
  ) %>%
  ungroup()

# Determine facet layout
tools <- sort(unique(df$tool))
modes <- unique(df$mode)
n_tools <- length(tools)
has_fair <- "fair" %in% modes

# Build violin + boxplot
p <- ggplot(df, aes(x = group_label, y = mod_ratio, fill = truth)) +
  geom_violin(scale = "width", alpha = 0.7, trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.5, linewidth = 0.3) +
  scale_fill_manual(values = c("positive" = "#E64B35", "negative" = "#4DBBD5")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    x = NULL,
    y = "Modification Ratio",
    fill = "Truth"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8)
  )

if (has_fair) {
  # Facet: tools as rows, mode as columns
  p <- p + facet_grid(tool ~ mode, scales = "free_x")
  fig_width <- 8
} else {
  # No fair mode (per-sample tools): facet by tool only
  p <- p + facet_wrap(~ tool, ncol = 1, scales = "free_x")
  fig_width <- 5
}

# Add comparison info to title if present
comparisons <- unique(df$comparison)
if (length(comparisons) == 1 && comparisons[1] != "unknown") {
  p <- p + ggtitle(paste0("Modification Ratio — ", comparisons[1]))
}

fig_height <- max(4, 2.5 * n_tools + 1)

ggsave(output_pdf, p, width = fig_width, height = fig_height, limitsize = FALSE)
message("Saved: ", output_pdf)
