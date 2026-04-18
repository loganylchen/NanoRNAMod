# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

NanoRNAMod is a Snakemake workflow for detecting RNA modifications (e.g., m6A) from Oxford Nanopore sequencing data. It orchestrates **comparison-based** detection tools that require both a "Native" sample (unmodified RNA) and a "Control" sample (IVT/synthetic RNA without modifications), then harmonizes their output into a unified format.

See `AGENTS.md` for code style guidelines and detailed command reference.

## Key Commands

```bash
# Lint and format
snakemake --lint --snakefile workflow/Snakefile
snakefmt workflow/Snakefile workflow/rules/*.smk

# Dry run (validate DAG without executing)
snakemake --use-conda --dry-run --cores 1

# Run full test suite
cd .test && snakemake --use-conda --use-singularity --show-failed-logs --cores 4 -k -p

# Run a specific rule or target
snakemake --use-conda -R <rule_name>
snakemake --use-conda --cores 2 -p "results/modifications/xpore/{comparison}/xpore_results.tsv"

# Generate HTML report
snakemake --report report.zip
```

## Architecture

### Data Flow

```
data/{sample}/fastq/pass.fq.gz        (raw ONT reads)
data/{sample}/blow5/nanopore.blow5    (raw signal)
        │
        ▼
{project}/results/fastq/{sample}.fq.gz    ← symlinked by link_fastq
{project}/results/blow5/{sample}.blow5    ← symlinked by link_blow5
        │
        ├─► minimap2 → BAM (genome + transcriptome)
        │         └─► f5c eventalign → collapse TSV (for baleen, nanocompore)
        │
        ├─► Modification detection tools (per comparison: {native}_{control})
        │         └─► {project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv
        │
        └─► QC, quantification, poly-A estimation, variant calling
```

### Sample Model

`samples.tsv` requires two conditions: `Native` (modified RNA) and `Control` (IVT/unmodified). Every pairwise combination of `native_sample × control_sample` forms a **comparison** wildcard (e.g., `A_F`). Most modification detection tools run per-comparison; tandemmod/directrm/m6atm/rnano run per-sample.

### Rule Organization

Rules are split by phase and prefixed accordingly:
- `prep_*` — data preparation (linking, alignment, filtering, indexing, signal alignment)
- `modetect_*` — one file per tool (xpore, nanocompore, baleen, pybaleen, differr, drummer, eligos2, epinano, tandemmod, directrm, m6atm, rnano, psipore, nanopsu, nanomud, penguin)
- `post_*` — output formatting and GTF annotation
- `qc_*` — QC modules (nanoplot, qualimap, samtools, nanocount, bcftools)
- `polya_*` — poly-A tail estimation (nanopolish)

### Tool Activation

Tools are enabled/disabled in `config/config.yaml` under `tools.<name>.activate`. Four tools (tandemmod, directrm, m6atm, rnano) are also commented out in the Snakefile `include` list — uncomment both the include and the config entry to activate them.

### Key Helper Functions (workflow/rules/common.smk)

- `get_container(tool_name)` — resolves container image, falling back to `default`
- `get_threads(tool_name, default)` — resolves thread count from config
- `get_final_output()` — builds the full list of target files based on activated tools
- `KEEP_OR_NOT(path)` — wraps path in `temp()` when `only_final_results: True`

### Output Path Convention

All results are rooted at `{project}/results/` (set via `config["project"]`). The project name acts as the top-level namespace, allowing multiple analyses in the same directory.

### Adding a New Tool

1. Create `workflow/envs/{tool}.yaml` (conda env)
2. Create `workflow/rules/modetect_{tool}.smk` with prep + detect rules
3. Create `workflow/scripts/{tool}_postprocess.py` to normalize output to TSV
4. Add `tools.{tool}.activate` to `config/config.yaml` and the schema
5. Add a `get_container("{tool}")` call in the rule and a thread entry in config
6. Add the tool's final output path to `get_final_output()` in `common.smk`
7. Uncomment or add the `include:` line in `workflow/Snakefile`

### Input Data Layout

```
data/
├── {sample}/
│   ├── fastq/pass.fq.gz  (or pass.fastq.gz)
│   └── blow5/nanopore.blow5  (or nanopore.drs.blow5)
└── ref.fa / ref.gtf      (reference files, paths set in config)
```
