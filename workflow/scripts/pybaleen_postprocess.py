#!/usr/bin/env python3
"""Post-process pybaleen site_results.tsv to standard NanoRNAMod format.

Standardizes column names:
  - contig -> transcript
  - position (unchanged)
  - pvalue -> p_value
  - padj -> adjusted_p_value
"""

import sys
import pandas as pd

input_file = snakemake.input.result
output_file = snakemake.output.result

df = pd.read_csv(input_file, sep='\t')

# Rename columns to standard format
df = df.rename(columns={
    'contig': 'transcript',
    'pvalue': 'p_value',
    'padj': 'adjusted_p_value',
})

# Write standardized output
df.to_csv(output_file, sep='\t', index=False)
