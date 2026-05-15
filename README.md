# NanoRNAMod

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.4.1-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/loganylchen/NanoRNAMod/workflows/Tests/badge.svg?branch=main)](https://github.com/loganylchen/NanoRNAMod/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for **RNA modification detection from Oxford Nanopore direct-RNA sequencing**. NanoRNAMod orchestrates comparison-based detection tools that require a paired **Native** sample (modified RNA from a biological source) and a **Control** sample (IVT/synthetic RNA without modifications), then harmonizes their outputs into a unified TSV format.

## Table of Contents

- [Features](#features)
- [Supported Tools](#supported-tools)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
  - [Sample Sheet](#sample-sheet)
  - [Reference Files](#reference-files)
  - [Tool Activation](#tool-activation)
  - [Thread Configuration](#thread-configuration)
  - [Container Configuration](#container-configuration)
  - [Tool Parameters](#tool-parameters)
- [Workflow Architecture](#workflow-architecture)
  - [Data Flow](#data-flow)
  - [Sample Model](#sample-model)
  - [Rule Organization](#rule-organization)
  - [Scheduling Priority](#scheduling-priority)
- [Output Files](#output-files)
- [Advanced Usage](#advanced-usage)
- [Adding a New Tool](#adding-a-new-tool)
- [Citation](#citation)

---

## Features

- **7 modification detection tools**: pybaleen (default tier-1), xpore, nanocompore, drummer, differr, eligos2, epinano
- **Unified output format**: per-tool results harmonized into `{tool}_results.tsv`
- **Containerized execution**: every rule pins a Docker/Singularity image for reproducibility
- **Depth-first scheduling**: priority ladder so finished comparisons release temp files before new ones start
- **Configurable tmp dir**: respect HPC scratch via `tmpdir:` config or `--default-resources tmpdir=...`
- **QC, quantification, poly-A, variant calling**: optional add-on lanes alongside modification detection

---

## Supported Tools

All currently shipped tools are **comparison-based** (require Native + Control). Per-sample tools were removed in favor of a single complete set — see `TODO.md` for the list of removed stubs and how to re-introduce them.

| Tool        | Method                                              | Container                                  |
|-------------|-----------------------------------------------------|--------------------------------------------|
| **pybaleen**| GPU-accelerated DTW + HMM (default-on)              | `btrspg/py-baleen-gpu:dev`                 |
| **xpore**   | Bayesian signal-intensity comparison                | `btrspg/xpore:2.1`                         |
| **nanocompore** | k-mer signal logit/GMM testing                  | `btrspg/nanocompore1:1.0.4`                |
| **drummer** | Base-call comparison with read-depth control        | `btrspg/drummer:92bb35a4…`                 |
| **differr** | Differential error-rate analysis                    | `btrspg/differr:0.2`                       |
| **eligos2** | Direct-RNA modification detection on aligned BAMs   | `btrspg/eligos2:2.1.0`                     |
| **epinano** | Mismatch + base-quality `sum_err` summary           | `btrspg/epinano:1.2.0`                     |

---

## Installation

### Prerequisites

- [Snakemake](https://snakemake.readthedocs.io/) ≥ 6.4.1 (8.x recommended)
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/) for environment management
- [Singularity / Apptainer](https://apptainer.org/) for `--use-singularity` (recommended)

### Clone the Repository

```bash
git clone https://github.com/loganylchen/NanoRNAMod.git
cd NanoRNAMod
```

### Install Snakemake

NanoRNAMod does not ship a single environment file; each rule pulls its own container or conda env. Install Snakemake itself with conda/mamba:

```bash
mamba create -n nanornamod -c conda-forge -c bioconda \
    "snakemake>=6.4.1" snakedeploy snakefmt
mamba activate nanornamod
```

If you prefer Singularity-only execution, you only need Snakemake; all tool dependencies come from the pinned Docker images listed under `containers:` in `config/config.yaml`.

---

## Quick Start

1. **Prepare your sample sheet** `config/samples.tsv` (TSV, at least `SampleName`, `Condition`):

   ```tsv
   SampleName	Condition
   A	Native
   B	Native
   F	Control
   ```

2. **Place input data** under `data/`:

   ```
   data/
   ├── A/
   │   ├── fastq/pass.fq.gz
   │   └── blow5/nanopore.blow5
   ├── B/
   │   ├── fastq/pass.fq.gz
   │   └── blow5/nanopore.blow5
   ├── F/
   │   ├── fastq/pass.fq.gz
   │   └── blow5/nanopore.blow5
   └── ref.fa
   ```

3. **Edit `config/config.yaml`** — set the project name, reference paths, and which tools to activate.

4. **Run the workflow**:

   ```bash
   # Validate the DAG without executing
   snakemake --use-singularity --dry-run --cores 1

   # Execute with Singularity containers (recommended)
   snakemake --use-singularity --cores 8 -p

   # Execute with conda envs (where containers are not pinned)
   snakemake --use-conda --cores 8 -p
   ```

---

## Configuration

### Sample Sheet

`config/samples.tsv` requires the following columns:

| Column        | Required | Description                                                        |
|---------------|----------|--------------------------------------------------------------------|
| `SampleName`  | yes      | Unique identifier; sample directory at `data/{SampleName}/`         |
| `Condition`   | yes      | `Native` (modified RNA) or `Control` (IVT/unmodified)               |

You must have at least one `Native` sample and at least one `Control` sample. Every `native × control` cross-product becomes a comparison wildcard (e.g. `A_F`, `B_F`).

### Reference Files

```yaml
reference:
  genome_fasta: ecoli_ref/ref.fa
  transcriptome_fasta: ecoli_ref/ref.fa   # may equal genome_fasta for prokaryotes
  transcriptome_gtf: ''                    # optional; leave '' to skip GTF annotation
```

### Tool Activation

Enable/disable in `config/config.yaml`:

```yaml
tools:
  xpore:       { activate: false }
  nanocompore: { activate: false }
  differr:     { activate: false }
  drummer:     { activate: false }
  eligos2:     { activate: false }
  epinano:     { activate: false }
  pybaleen:    { activate: true }   # default-on
```

Only the seven tools above are recognized; setting `activate: true` on any other name is silently ignored.

### Thread Configuration

```yaml
threads:
  default: 1
  minimap2: 4
  f5c: 4
  slow5tools: 4
  xpore: 4
  nanocompore: 4
  differr: 4
  drummer: 4
  eligos2: 4
  epinano: 4
  pybaleen: 4
```

`get_threads(tool, default)` in `common.smk` resolves each rule's thread count from this map, falling back to `threads.default` then to the rule's hard-coded default.

### Container Configuration

```yaml
containers:
  default: "docker://condaforge/mambaforge:22.11.1-4"
  minimap2: "docker://btrspg/minimap2:2.28"
  f5c: "docker://btrspg/f5c:1.5"
  nanocompore: "docker://btrspg/nanocompore1:1.0.4"
  xpore: "docker://btrspg/xpore:2.1"
  differr: "docker://btrspg/differr:0.2"
  eligos2: "docker://btrspg/eligos2:2.1.0"
  drummer: "docker://btrspg/drummer:92bb35a4a2b22ff304f5e4bcbc9fa6985f18a12e"
  epinano: "docker://btrspg/epinano:1.2.0"
  pybaleen: "docker://btrspg/py-baleen-gpu:dev"
  # ... see config/config.yaml for the full list
```

`get_container(tool)` returns the matching image; falls back to `containers.default` when an entry is empty or missing.

### Tool Parameters

Per-tool CLI flags live under `params:` in `config/config.yaml`:

```yaml
params:
  minimap2_transcriptome: " -ax map-ont -L --secondary=no -N 10 -p 0 "
  minimap2_genome:        " -ax splice -u f -k 14 -G 500000 --secondary=no "
  samtools_filtering:     " -F 2324 -q 10 --min-qlen 80 "
  f5c_eventalign_full:    "--print-read-names --scale-events --samples --rna "
  f5c_eventalign_simple:  "--signal-index --scale-events --rna "
  nanocompore: "--logit --overwrite "
  xpore: ""
  differr: " -f 2  --median-expr-threshold 0  --min-expr-threshold 0 "
  eligos2: "--max_depth 5000000 --esb 0 --oddR 0 --pval 1"
  epinano: '-f sum_err'
  drummer: ' -z 0.1 -p 2 '
  pybaleen: ""
```

---

## Workflow Architecture

### Data Flow

```
data/{sample}/fastq/pass.fq.gz        (raw ONT reads)
data/{sample}/blow5/nanopore.blow5    (raw signal, slow5 format)
        │
        ▼
{project}/results/fastq/{sample}.fq.gz       ← link_fastq
{project}/results/blow5/{sample}.blow5       ← link_blow5
        │
        ├─► minimap2 → BAM (genome + transcriptome)
        │         └─► f5c eventalign → collapse TSV (xpore, nanocompore, pybaleen)
        │
        ├─► Modification detection per comparison {native}_{control}
        │         └─► {project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv
        │
        └─► Optional: QC, NanoCount quantification, nanopolish polya, bcftools variants
```

### Sample Model

The `Condition` column drives the comparison wildcard. With Native = {A, B} and Control = {F}, the workflow produces results for the comparisons `A_F` and `B_F`. Each comparison runs every activated tool independently.

### Rule Organization

| Prefix       | Purpose                                            |
|--------------|----------------------------------------------------|
| `prep_*`     | Linking, alignment (minimap2), filtering, indexing, f5c eventalign |
| `modetect_*` | One file per tool: pybaleen, xpore, nanocompore, differr, drummer, eligos2, epinano |
| `post_*`     | Normalize each tool's output into the unified TSV  |
| `qc_*`       | NanoPlot, Qualimap, samtools stats, NanoCount, bcftools |
| `polya_*`    | nanopolish polya estimation                        |

### Scheduling Priority

Rules use `priority` to enforce **depth-first** scheduling — finishing a comparison end-to-end frees its `temp()` intermediates before another comparison starts. Higher number = scheduled first.

| Tier | Range  | Tools / Stage                                                          |
|------|--------|------------------------------------------------------------------------|
| 1    | 96–100 | pybaleen and its dependencies (link_fastq=96, minimap2=97, f5c_index=98, eventalign=89/82, pybaleen_run=99, post_pybaleen=100) |
| 2    | 89–95  | xpore                                                                  |
| 3    | 82–86  | nanocompore                                                            |
| 4    | 76–78  | drummer                                                                |
| 5    | 74–75  | differr                                                                |
| 6    | 70–72  | eligos2                                                                |
| 7    | 66–68  | epinano                                                                |

This keeps the disk footprint bounded even when many comparisons run concurrently.

---

## Output Files

```
{project}/results/
├── workflow_version.json
├── alignments/{sample}.bam{,.bai}
├── alignments/{sample}_filtered.bam{,.bai}
├── quantification/{sample}.tx_counts.tsv
├── modifications/{tool}/{native}_{control}/{tool}_results.tsv   # unified output
└── qc/, polya/, variants/                                       # opt-in lanes
```

Each `{tool}_results.tsv` is a TSV produced by `workflow/scripts/format.py` (or per-tool `*_postprocess.py` for pybaleen). Empty-but-successful runs produce zero-byte sentinel files rather than failing the DAG.

---

## Advanced Usage

```bash
# Lint
snakemake --lint --snakefile workflow/Snakefile
snakefmt workflow/Snakefile workflow/rules/*.smk

# Force re-run a single rule
snakemake --use-singularity -R post_pybaleen --cores 4

# Build only one specific output
snakemake --use-singularity --cores 4 \
    "ecoliRNA/results/modifications/pybaleen/A_F/pybaleen_results.tsv"

# Generate the HTML report
snakemake --report report.zip

# Run the smoke test
cd .test && snakemake --use-conda --use-singularity \
    --show-failed-logs --cores 4 -k -p
```

### Tmp Directory

Many HPC nodes have small `/tmp`. NanoRNAMod sets `TMPDIR` for all subprocesses and registers a Snakemake default resource. Override per run via:

```bash
snakemake --default-resources "tmpdir='/scratch/$USER/tmp'" ...
```

or persistently via `tmpdir:` in `config/config.yaml`.

---

## Adding a New Tool

1. `workflow/envs/{tool}.yaml` — conda env for the tool (if not pulling a container)
2. `workflow/rules/modetect_{tool}.smk` — prep + detect rules, set `priority:` according to the tier scheme above
3. `workflow/scripts/{tool}_postprocess.py` — normalize the tool's output into the unified TSV
4. `config/config.yaml` — add entries under `params:`, `threads:`, `containers:`, `tools:`
5. `workflow/rules/common.smk` — append the tool name to `PER_COMPARISON_TOOLS` or `PER_SAMPLE_TOOLS`
6. `workflow/rules/post_format.smk` — add a `post_{tool}` rule (or a branch in `format.py`)
7. `workflow/Snakefile` — add an `include: "rules/modetect_{tool}.smk"` line

See `TODO.md` for the list of tools previously removed and the same recipe applied in reverse.

---

## Citation

If you use this workflow in a paper, please cite:

- This repository: https://github.com/loganylchen/NanoRNAMod
- Snakemake: Mölder et al., 2021, [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2)
- The specific modification-detection tools you activate (xpore, nanocompore, etc.)

---

## License

MIT — see `LICENSE`.

---

## Acknowledgments

This workflow integrates and harmonizes the output of multiple independent RNA modification detection tools. Thanks to the developers of each underlying tool for making their work available.
