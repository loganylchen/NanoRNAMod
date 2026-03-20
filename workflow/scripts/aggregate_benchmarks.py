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
# Aggregate benchmark files into a single summary table
# ─────────────────────────────────────────────────────────────────────────────

rows = []

for f in sorted(benchmark_files):
    try:
        df = pd.read_csv(f, sep='\t')
        # Extract tool and sample from filename
        stem = os.path.basename(f).replace(".benchmark.txt", "")
        tool = stem.rsplit(".", 1)[-1]
        sample = stem.rsplit(".", 1)[0]
        df['tool'] = tool
        df['sample'] = sample
        rows.append(df)
    except Exception as e:
        print(f"Warning: could not parse {f}: {e}")

if rows:
    df = pd.concat(rows, ignore_index=True)
    df = df.sort_values(["tool", "sample"])
else:
    df = pd.DataFrame()

df.to_csv(output_file, sep='\t', index=False)
