# Snakemake workflow: `NanoRNAMod`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/loganylchen/NanoRNAMod/workflows/Tests/badge.svg?branch=main)](https://github.com/loganylchen/NanoRNAMod/actions?
query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for `RNA modification detection from Nanopore sequencing data`. Especically for comparison based methods.


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=loganylchen%2FNanoRNAMod).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <repo>sitory and its DOI (see above).

## Container and Thread Configuration

### Container Configuration | 容器配置

All rules in this workflow can use custom container images. Configure them in `config/config.yaml`:

所有规则都可以使用自定义容器镜像。在 `config/config.yaml` 中配置：

```yaml
containers:
  default: "docker://condaforge/mambaforge:22.11.1-4"  # Default container for all tools
  xpore: "docker://myregistry/xpore:latest"            # Custom container for xpore
  minimap2: "docker://biocontainers/minimap2:2.24"     # Custom container for minimap2
  nanocompore: ""                                       # Empty string uses default container
```

**How it works | 工作原理:**
- Each rule calls `get_container("tool_name")` to get its container image
- If a tool's container is empty `""` or not specified, it uses the `default` container
- Supports Docker (`docker://`) and Singularity Hub (`shub://`) formats

**Running with containers | 使用容器运行:**
```bash
snakemake --use-conda --use-singularity --cores 4
```

### Thread Configuration | 线程配置

Configure thread counts for each tool in `config/config.yaml`:

在 `config/config.yaml` 中配置每个工具的线程数：

```yaml
threads:
  default: 1      # Default threads when not specified
  minimap2: 8     # Use 8 threads for minimap2
  xpore: 16       # Use 16 threads for xpore
  nanocompore: 12 # Use 12 threads for nanocompore
```

**How it works | 工作原理:**
- Each rule calls `get_threads("tool_name", default_value)` to get thread count
- If not specified in config, uses the `default` value
- If `default` is not set, uses the rule's default value

### Available Tools | 可用工具

The following tools support custom container and thread configuration:

以下工具支持自定义容器和线程配置：

- **Alignment**: minimap2
- **Signal processing**: f5c, slow5tools
- **Modification detection**: xpore, nanocompore, baleen, differr, drummer, eligos2, epinano, tandemmod, directrm, m6atm, rnano
- **QC tools**: samtools, nanoplot, nanocount, qualimap, bcftools, nanopolish
