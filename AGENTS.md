# AGENTS.md

This file contains guidelines for agentic coding agents working on the NanoRNAMod repository.

## Build/Lint/Test Commands

### Snakemake Workflow Commands
```bash
snakemake --lint --snakefile workflow/Snakefile
snakefmt workflow/Snakefile workflow/rules/*.smk
cd .test && snakemake --use-conda --use-singularity --show-failed-logs --cores 4 -k -p
snakemake --use-conda -R <rule_name>
snakemake --use-conda --dry-run --cores 1
snakemake --report report.zip
```

### Running Individual Tests
```bash
cd .test && snakemake --use-conda --cores 2 -p "results/modifications/xpore/{comparison}/xpore_results.tsv"
snakemake --use-conda --cores 2 -R xpore_dataprep
snakemake --use-conda -F <target_file>
```

### GitHub Actions
- **Formatting**: `github/super-linter@v4` with snakefmt
- **Linting**: `snakemake/snakemake-github-action` with `--lint`
- **Testing**: Full workflow with conda/singularity

## Code Style Guidelines

### Python Scripts (workflow/scripts/)

**Imports**: Standard library first, then third-party
```python
import os
import sys
import pandas as pd
import pysam
from snakemake.shell import shell
```

**Snakemake Integration**: Access via `snakemake.input[0]`, `snakemake.output[0]`, `snakemake.params.param_name`, `snakemake.threads`. Redirect stdout/stderr:
```python
sys.stdout = open(snakemake.log[0], "w")
sys.stderr = open(snakemake.log[1], "w")
```

**File Paths**: Use `os.path.dirname()`, `os.path.abspath()`, `os.makedirs(path, exist_ok=True)`, f-strings: `f"{outdir}/{transcript}.bam"`

**Shell Commands**: Use `subprocess.Popen()` or `snakemake.shell()`:
```python
subprocess.Popen(f"command", shell=True, stderr=sys.stdout, stdout=sys.stdout).communicate()
shell("command {input} {output} 2>{log}")
```

**Data Processing**: pandas (`df = pd.read_csv(file, sep='\t')`), `df.sort_values()`, `df.to_csv()`, `collections.defaultdict`

**Multiprocessing**: `multiprocessing.Pool()` for parallel tasks

**No Type Hints**: This project does not use Python type hints

### Snakemake Rules (workflow/rules/)

**Rule Naming**
- Use lowercase_with_underscores: `rule xpore_dataprep:`
- Use descriptive names that indicate function

**Input/Output/Log**
- Use consistent naming: `{project}/results/{tool}/{sample}_file.ext`
- Use `temp()` for intermediate files that should be deleted: `temp("{path}")`
- Always define log files: `log: "logs/{project}/{rule}/{sample}.log"`
- Use `benchmark:` for performance tracking

**Resources**
- Define memory: `resources: mem_mb = 1024 * 50`
- Set threads from config: `threads: config["threads"]["tool"]`
- Set priority: `priority: 10`

**Conda Environments**
- Reference env files: `conda: "../envs/tool.yaml"`

**Shell Commands**
- Use parameter expansion: `{input}`, `{output}`, `{threads}`, `{log}`
- Redirect stderr to log: `2>{log}`
- Chain commands with `&&` for dependency management
- Use spaces around operators in shell commands

**Wildcards**
- Define in wildcard_constraints at top of Snakefile
- Use patterns: `sample="[\da-zA-Z]+_?[0-9]*"`

**Script Sections**
- Use Python scripts for complex logic instead of shell: `script: "../scripts/tool.py"`

### Directory Structure
```
workflow/
├── Snakefile              # Main workflow with includes
├── rules/                 # Individual rule files (.smk)
├── scripts/               # Python helper scripts
└── envs/                  # Conda environment YAML files
config/
├── config.yaml            # Main configuration
└── samples.tsv            # Sample definitions
.test/
├── config/                # Test configurations
└── data/                  # Test data
```

### General Guidelines

**Error Handling**
- Use sys.exit(1) for critical errors in common.smk
- Log errors to stderr for shell commands
- Validate input paths exist before use

**Configuration**
- Use `config["key"]` for accessing configuration values
- Validate config with schemas in `workflow/schemas/`
- Define defaults in config files

**Logging**
- All rules should have log files
- Python scripts redirect stdout/stderr to log files
- Print progress messages for long-running operations

**Testing**
- Place test data in `.test/data/`
- Test configs in `.test/config/`
- Use sample data that fits reasonable runtime for CI

**No Comments**: Follow project convention of minimal inline comments unless necessary for clarity
