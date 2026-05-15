#!/usr/bin/env python3
"""Post-process pybaleen site_results.tsv to standard NanoRNAMod format.

Standardizes column names:
  - contig -> transcript
  - position (unchanged)
  - pvalue -> p_value
  - padj -> adjusted_p_value
"""

import os
import sys
import pandas as pd

input_file = snakemake.input.result
output_file = snakemake.output.result

os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)

# Empty-input safety net: pybaleen ran successfully but produced no sites
# (e.g. depth/quality filter rejected everything). Propagate as empty output.
if not (os.path.isfile(input_file) and os.path.getsize(input_file) > 0):
    open(output_file, "w").close()
    sys.exit(0)

df = pd.read_csv(input_file, sep='\t')

# Rename columns to standard format
df = df.rename(columns={
    'contig': 'transcript',
    'pvalue': 'p_value',
    'padj': 'adjusted_p_value',
})

# Write standardized output
df.to_csv(output_file, sep='\t', index=False)
