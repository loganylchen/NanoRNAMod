"""
Aggregate Snakemake benchmark files into a single resource summary TSV.

Snakemake benchmark file columns (tab-separated, first row is header):
    s           - wall-clock time in seconds
    h:m:s       - wall-clock time as HH:MM:SS
    max_rss     - maximum resident set size (MB)
    max_vms     - maximum virtual memory size (MB)
    max_uss     - maximum unique set size (MB)
    max_pss     - maximum proportional set size (MB)
    io_in       - total I/O read (MB)
    io_out      - total I/O write (MB)
    mean_load   - mean CPU load (%)
    cpu_time    - total CPU time (seconds)

Filename convention:  {stem}.benchmark.txt
    where stem is typically  {sample}.{tool}  or  {native}_{control}.{tool}
    e.g. SampleA.tandemmod.benchmark.txt  →  sample="SampleA", tool="tandemmod"
         SampleA_SampleB.eligos2.benchmark.txt  →  sample="SampleA_SampleB", tool="eligos2"
"""

import os
import glob
import pandas as pd

# Snakemake passes input/output/params via the snakemake object.
benchmark_files = snakemake.input          # list of .benchmark.txt paths
output_file = snakemake.output[0]
benchmark_dir = snakemake.params.benchmark_dir

# ── If input is empty (no benchmarks exist yet), fall back to a glob ──────────
if not benchmark_files:
    benchmark_files = glob.glob(os.path.join(benchmark_dir, "*.benchmark.txt"))

# ─────────────────────────────────────────────────────────────────────────────
# TODO: Implement the aggregation logic below (~10 lines).
#
# Steps to complete:
#   1. Parse each benchmark file with pd.read_csv(f, sep='\t')
#   2. Extract `tool` and `sample` from the filename stem:
#        stem = os.path.basename(f).replace(".benchmark.txt", "")
#        tool  = stem.rsplit(".", 1)[-1]      # last dot-separated token
#        sample = stem.rsplit(".", 1)[0]      # everything before the last dot
#   3. Add "tool" and "sample" columns to each DataFrame
#   4. Concatenate all DataFrames with pd.concat(...)
#   5. Sort by ["tool", "sample"] and write to output_file with sep='\t', index=False
#
# Trade-offs to consider:
#   - rsplit(".", 1) handles tool names with no dots; stems like
#     "SampleA_SampleB.eligos2" correctly yield sample="SampleA_SampleB"
#   - You may want to add error handling for malformed or empty files
# ─────────────────────────────────────────────────────────────────────────────

rows = []

for f in sorted(benchmark_files):
    # TODO: parse the file and extract tool/sample, then append to rows
    pass

# TODO: build a DataFrame from rows and write to output_file
# Example final lines:
#   df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
#   df = df.sort_values(["tool", "sample"]) if not df.empty else df
#   df.to_csv(output_file, sep='\t', index=False)
