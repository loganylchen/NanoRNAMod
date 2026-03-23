# NanoRNAMod

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/loganylchen/NanoRNAMod/workflows/Tests/badge.svg?branch=main)](https://github.com/loganylchen/NanoRNAMod/actions?query=branch%3Amain+workflow%3ATests)

A comprehensive Snakemake workflow for **RNA modification detection from Oxford Nanopore sequencing data**. This workflow orchestrates comparison-based detection tools that require both a "Native" sample (modified RNA from biological samples) and a "Control" sample (IVT/synthetic RNA without modifications), then harmonizes their output into a unified format with integrated benchmarking capabilities.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
  - [Sample Sheet](#sample-sheet)
  - [Reference Files](#reference-files)
  - [Tool Activation](#tool-activation)
  - [Thread Configuration](#thread-configuration)
  - [Container Configuration](#container-configuration)
  - [Benchmarking Configuration](#benchmarking-configuration)
- [Workflow Architecture](#workflow-architecture)
  - [Data Flow](#data-flow)
  - [Rule Organization](#rule-organization)
- [Supported Tools](#supported-tools)
  - [Comparison-Based Tools](#comparison-based-tools)
  - [Per-Sample Tools](#per-sample-tools)
- [Output Files](#output-files)
  - [Modification Results](#modification-results)
  - [Benchmarking Results](#benchmarking-results)
- [Benchmarking Module](#benchmarking-module)
  - [P-value Transformation](#p-value-transformation)
  - [Evaluation Metrics](#evaluation-metrics)
  - [Multi-Threshold Analysis](#multi-threshold-analysis)
- [Advanced Usage](#advanced-usage)
- [Citation](#citation)

---

## Features

- **16+ Modification Detection Tools**: Integrates xpore, nanocompore, baleen, pybaleen, differr, drummer, eligos2, epinano, tandemmod, directrm, m6atm, rnano, psipore, nanopsu, nanomud, and penguin
- **Unified Output Format**: All tool outputs are harmonized into a standardized TSV format
- **Containerized Execution**: Full Docker/Singularity support for reproducibility
- **Flexible Configuration**: Per-tool thread counts, containers, and parameters
- **Comprehensive Benchmarking**: Built-in accuracy evaluation with precision, recall, F1, AUPRC, AUROC, and more
- **P-value Transformation**: Automatic -log10 transformation for uniform score comparison
- **Interactive Reports**: HTML benchmarking reports with visualization plots

---

## Installation

### Prerequisites

- [Snakemake](https://snakemake.readthedocs.io/) >= 6.3.0
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)
- [Singularity](https://sylabs.io/singularity/) (optional, for containerized execution)

### Clone the Repository

```bash
git clone https://github.com/loganylchen/NanoRNAMod.git
cd NanoRNAMod
```

### Install Dependencies

```bash
# Using conda
conda env create -f workflow/envs/main.yaml
conda activate nanornamod

# Or using mamba (faster)
mamba env create -f workflow/envs/main.yaml
mamba activate nanornamod
```

---

## Quick Start

1. **Prepare your sample sheet** (`config/samples.tsv`):
   ```tsv
   SampleName	Condition	Directory
   SampleA	Native	path/to/sampleA
   SampleB	Native	path/to/sampleB
   SampleC	Control	path/to/sampleC
   ```

2. **Configure reference files** in `config/config.yaml`:
   ```yaml
   reference:
     genome_fasta: data/genome.fa
     transcriptome_fasta: data/transcriptome.fa
     transcriptome_gtf: data/annotations.gtf
   ```

3. **Run the workflow**:
   ```bash
   # Dry run to check configuration
   snakemake --use-conda --dry-run --cores 1

   # Run with conda environments
   snakemake --use-conda --cores 4

   # Run with singularity containers (recommended for reproducibility)
   snakemake --use-singularity --cores 4
   ```

---

## Configuration

### Sample Sheet

The sample sheet (`config/samples.tsv`) defines your input samples. It must contain at least three columns:

| Column | Description |
|--------|-------------|
| `SampleName` | Unique identifier for each sample |
| `Condition` | Either `Native` (modified RNA) or `Control` (IVT/unmodified) |
| `Directory` | Path to the sample directory containing FASTQ/BLOW5 files |

**Example:**
```tsv
SampleName	Condition	Directory
HeLa_Rep1	Native	data/HeLa_Rep1
HeLa_Rep2	Native	data/HeLa_Rep2
IVT_Rep1	Control	data/IVT_Rep1
IVT_Rep2	Control	data/IVT_Rep2
```

**Input Data Layout:**
```
data/
├── {sample}/
│   ├── fastq/pass.fq.gz      (or pass.fastq.gz)
│   └── blow5/nanopore.blow5  (or nanopore.drs.blow5)
└── ref.fa / ref.gtf          (reference files)
```

### Reference Files

Configure your reference files in `config/config.yaml`:

```yaml
reference:
  genome_fasta: data/ref.fa           # Genome reference
  transcriptome_fasta: data/ref.fa    # Transcriptome reference (can be same as genome)
  transcriptome_gtf: data/ref.gtf     # Gene annotations in GTF format
```

### Tool Activation

Enable or disable tools in `config/config.yaml`:

```yaml
tools:
  xpore:
    activate: True
  nanocompore:
    activate: True
  baleen:
    activate: True
  # ... other tools
```

**Note:** Four tools (tandemmod, directrm, m6atm, rnano) may also need to be uncommented in the Snakefile `include` list to activate them.

### Thread Configuration

Configure thread counts for each tool to optimize parallelization:

```yaml
threads:
  default: 1           # Default threads when not specified
  minimap2: 8          # Alignment benefits from more threads
  nanocompore: 4       # Moderate parallelization
  xpore: 4
  baleen: 4
  # ... other tools
```

**How it works:**
- Each rule calls `get_threads("tool_name", default_value)` to get thread count
- If not specified in config, uses the provided default value
- Adjust based on your system's CPU cores

### Container Configuration

Configure container images for reproducibility:

```yaml
containers:
  default: "docker://condaforge/mambaforge:22.11.1-4"
  minimap2: "docker://btrspg/minimap2:2.28"
  nanocompore: "docker://btrspg/nanocompore1:1.0.4"
  # ... other tools
```

**How it works:**
- Each rule calls `get_container("tool_name")` to get its container image
- If a tool's container is empty `""` or not specified, it uses the `default` container
- Supports Docker (`docker://`) and Singularity Hub (`shub://`) formats

**Running with containers:**
```bash
snakemake --use-conda --use-singularity --cores 4
```

### Benchmarking Configuration

Configure accuracy benchmarking against a ground truth set:

```yaml
benchmark:
  truth_set: "config/benchmark_truth.tsv"  # Path to ground truth TSV
  window: [0, 1, 2, 5]                     # Position tolerance windows (nt)
  n_thresholds: 50                          # Number of thresholds to evaluate
  custom_thresholds: []                     # Optional custom thresholds
```

**Truth Set Format:**
```tsv
transcript	position	modification_type	label
ENST001	        1234	        m6A	        positive
ENST001	        5678	        m6A	        -
ENST002	        100	        m6A	                # No label = positive (backward compatible)
```

**Label Options:**
- `positive` - Known modification site
- `-` - Known negative site (unmodified)

---

## Workflow Architecture

### Data Flow

```
data/{sample}/fastq/pass.fq.gz        (raw ONT reads)
data/{sample}/blow5/nanopore.blow5    (raw signal)
        │
        ▼
{project}/results/fastq/{sample}.fq.gz    ← symlinked by link_fastq
{project}/results/blow5/{sample}.blow5    ← symlinked by link_blow5
        │
        ├─► minimap2 → BAM (genome + transcriptome alignment)
        │         └─► f5c eventalign → collapse TSV (for signal-based tools)
        │
        ├─► Modification detection tools (per comparison: {native}_{control})
        │         └─► {project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv
        │
        └─► QC, quantification, poly-A estimation, variant calling
```

### Sample Model

`samples.tsv` requires two conditions:
- **Native**: Modified RNA from biological samples
- **Control**: IVT (in vitro transcribed) or synthetic RNA without modifications

Every pairwise combination of `native_sample × control_sample` forms a **comparison** wildcard (e.g., `A_F`, `B_F`). Most modification detection tools run per-comparison.

### Rule Organization

Rules are split by phase and prefixed accordingly:

| Prefix | Purpose | Examples |
|--------|---------|----------|
| `prep_*` | Data preparation | linking, alignment, filtering, indexing, signal alignment |
| `modetect_*` | Modification detection | xpore, nanocompore, baleen, etc. (one file per tool) |
| `post_*` | Output formatting | GTF annotation, result aggregation |
| `qc_*` | Quality control | nanoplot, qualimap, samtools, nanocount, bcftools |
| `polya_*` | Poly-A tail estimation | nanopolish |
| `benchmark_*` | Accuracy evaluation | multi-threshold analysis, visualization |

---

## Supported Tools

### Comparison-Based Tools

These tools require paired Native + Control samples:

| Tool | Description | Method |
|------|-------------|--------|
| **xpore** | Differential RNA modification detection | Signal intensity comparison |
| **nanocompore** | Nanopore RNA modification calling | k-mer level signal analysis |
| **baleen** | Modification detection using neural networks | Signal alignment + ML |
| **pybaleen** | Python reimplementation of baleen | GPU-accelerated |
| **differr** | Differential error rate analysis | Base-calling error patterns |
| **drummer** | Detection of modified bases | Signal comparison |
| **eligos2** | Epitranscriptional modification detection | RNA-specific features |
| **epinano** | Epitranscriptome analysis | Mismatch + signal features |
| **psipore** | Pseudouridine detection | Signal-based |

### Per-Sample Tools

These tools run on individual samples without requiring a control:

| Tool | Description | Target |
|------|-------------|--------|
| **tandemmod** | Tandem modification detection | Multiple modifications |
| **directrm** | Direct RNA modification detection | Multiple modifications |
| **m6atm** | m6A detection tool | m6A |
| **rnano** | RNA modification detection | Multiple modifications |
| **nanopsu** | Pseudouridine detection | Psi |
| **nanomud** | Modification detection | Psi, m1Psi |
| **penguin** | Modification detection | Multiple |

### Tool Parameters

Configure tool-specific parameters in `config/config.yaml`:

```yaml
params:
  minimap2_transcriptome: " -ax map-ont -L --secondary=no -N 10 -p 0 "
  minimap2_genome: " -ax splice -u f -k 14 -G 500000 --secondary=no "
  nanocompore: "--logit --overwrite"
  xpore: ""
  baleen_modcall: ""
  differr: " -f 2 --median-expr-threshold 0 --min-expr-threshold 0 "
  eligos2: "--max_depth 5000000 --esb 0 --oddR 0 --pval 1"
  epinano: " -t 1 -c 5 -f sum_err -d 0.01 "
  drummer: " -z 0.1 -p 1 "
  tandemmod: "--model multi --threshold 0.5"
  directrm: "--mods all --min-prob 0.7"
  m6atm: "--stoichiometry --threshold 0.6"
  pybaleen: "--padding 0 --min-depth 15"
  # ... other tools
```

---

## Output Files

### Modification Results

All modification detection results are stored in:
```
{project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv
```

**Unified Output Format:**
| Column | Description |
|--------|-------------|
| `transcript` | Transcript ID |
| `position` | Genomic position (1-based) |
| `ref_kmer` | Reference k-mer sequence |
| `coverage` | Read coverage at this position |
| `score` | Tool-specific significance score |
| `mod_type` | Predicted modification type |
| `comparison` | Native vs Control comparison ID |

### Benchmarking Results

When benchmarking is enabled, results are stored in:
```
{project}/results/benchmarks/
├── accuracy_summary.tsv             # Per-modification metrics
├── accuracy_summary_overall.tsv     # Aggregated metrics across all tools
├── optimal_thresholds.tsv           # Best score thresholds per tool
├── threshold_evaluation.tsv         # Metrics at each evaluated threshold
├── optimal_thresholds_detailed.tsv  # Comprehensive optimal threshold analysis
├── score_distributions.tsv          # Score statistics per tool
├── detailed_predictions.tsv         # Site-by-site prediction analysis
├── detailed_truth.tsv               # Site-by-site truth coverage
└── viz/
    └── benchmark_report.html        # Interactive HTML report with plots
```

---

## Benchmarking Module

### P-value Transformation

The benchmarking module automatically detects and transforms p-value columns to ensure uniform "higher = better" score comparison:

```python
# Automatic detection of p-value columns
def is_pvalue_column(col_name):
    """Detect if a column contains p-values based on name patterns."""
    patterns = ['pval', 'p_value', 'pvalue', 'p.value', 'adjusted_pval', 'padj', 'fdr']
    return any(p in col_name.lower() for p in patterns)

# -log10 transformation for p-values
def transform_pvalue_to_log10(scores):
    """Transform p-values to -log10(p-values) for uniform comparison."""
    return -np.log10(scores.clip(lower=1e-300))
```

**Why this matters:**
- Raw p-values: lower = more significant (requires reverse sorting)
- Transformed (-log10): higher = more significant (consistent with other scores)
- This ensures all tools are evaluated consistently regardless of their output format

### Evaluation Metrics

The benchmarking module computes comprehensive metrics at multiple score thresholds:

| Metric | Description | Formula |
|--------|-------------|---------|
| **Precision** | True positives / predicted positives | TP / (TP + FP) |
| **Recall** | True positives / actual positives | TP / (TP + FN) |
| **F1 Score** | Harmonic mean of precision and recall | 2 × (P × R) / (P + R) |
| **Specificity** | True negatives / actual negatives | TN / (TN + FP) |
| **MCC** | Matthews Correlation Coefficient | Balanced measure for imbalanced data |
| **AUPRC** | Area Under Precision-Recall Curve | Threshold-independent |
| **AUROC** | Area Under ROC Curve | Threshold-independent |

### Multi-Threshold Analysis

The workflow evaluates performance across multiple thresholds:

1. **Threshold Generation**: Creates `n_thresholds` evenly spaced thresholds from min to max scores
2. **Per-Window Evaluation**: Computes metrics at each position tolerance window (0, 1, 2, 5 nt)
3. **Optimal Threshold Selection**: Identifies threshold maximizing F1 score
4. **Score Distribution Analysis**: Reports score statistics for positive/negative sites

**Output Columns:**
| Column | Description |
|--------|-------------|
| `tool` | Tool name |
| `comparison` | Sample comparison ID |
| `modification_type` | Type of modification (e.g., m6A, Psi) |
| `window` | Position tolerance in nucleotides |
| `threshold` | Score threshold evaluated |
| `precision` | Precision at this threshold |
| `recall` | Recall at this threshold |
| `f1` | F1 score at this threshold |
| `tp/fp/fn/tn` | Confusion matrix counts |
| `auprc/auroc` | Area under curves |
| `original_score_column` | Which score column was used |
| `original_threshold` | Threshold in original (pre-transformed) scale |

---

## Advanced Usage

### Run Specific Rules

```bash
# Run a specific rule
snakemake --use-conda -R <rule_name>

# Generate a specific output file
snakemake --use-conda --cores 2 "results/modifications/xpore/A_F/xpore_results.tsv"
```

### Lint and Format

```bash
# Lint the workflow
snakemake --lint --snakefile workflow/Snakefile

# Format workflow files
snakefmt workflow/Snakefile workflow/rules/*.smk
```

### Generate HTML Report

```bash
# Create a comprehensive HTML report
snakemake --report report.zip
```

### Dry Run

```bash
# Validate DAG without executing
snakemake --use-conda --dry-run --cores 1
```

### Run Tests

```bash
# Run the test suite
cd .test && snakemake --use-conda --use-singularity --show-failed-logs --cores 4 -k -p
```

### Adding a New Tool

1. Create `workflow/envs/{tool}.yaml` (conda environment)
2. Create `workflow/rules/modetect_{tool}.smk` with prep + detect rules
3. Create `workflow/scripts/{tool}_postprocess.py` to normalize output to TSV
4. Add `tools.{tool}.activate` to `config/config.yaml` and the schema
5. Add a `get_container("{tool}")` call in the rule and a thread entry in config
6. Add the tool's final output path to `get_final_output()` in `common.smk`
7. Uncomment or add the `include:` line in `workflow/Snakefile`

---

## Citation

If you use this workflow in a paper, please cite:

- The URL of this repository: https://github.com/loganylchen/NanoRNAMod
- The DOI of the repository
- Snakemake: Mölder et al., 2021, [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2)

Also cite the specific modification detection tools you use in your analysis.

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

---

## Acknowledgments

This workflow integrates and harmonizes the output of multiple RNA modification detection tools. We thank all the developers of the individual tools for their contributions to the field.
