# Implementation Plan: Adding Missing Nanopore RNA Modification Tools

## Executive Summary

This document outlines the systematic implementation of **10 high-priority tools** identified in the literature review to expand the NanoRNAMod workflow from 7 to 17 total tools. The implementation is organized into **4 phases** based on complexity, dependencies, and priority.

**Timeline**: 8-12 weeks total (2-3 weeks per phase)
**Total Tools to Add**: 10 tools
**Estimated Lines of Code**: ~3,000-4,000
**New Environments**: 10 conda env files
**New Rule Files**: 10 .smk files
**New Scripts**: 15-20 helper scripts

---

## Table of Contents

1. [Overview and Architecture](#overview-and-architecture)
2. [Phase-by-Phase Implementation Plan](#phase-by-phase-implementation-plan)
3. [Tool-Specific Implementation Details](#tool-specific-implementation-details)
4. [Technical Specifications](#technical-specifications)
5. [Testing and Validation Strategy](#testing-and-validation-strategy)
6. [Dependency Management](#dependency-management)
7. [Configuration Updates](#configuration-updates)
8. [Documentation Requirements](#documentation-requirements)
9. [Timeline and Milestones](#timeline-and-milestones)
10. [Risk Assessment and Mitigation](#risk-assessment-and-mitigation)

---

## Overview and Architecture

### Current Architecture

```
NanoRNAMod Workflow
├── Data Preparation (minimap2, filtering, indexing)
├── Eventalign (f5c) - Required for: baleen, nanocompore
├── Modification Detection (7 tools)
│   ├── xpore (comparative)
│   ├── baleen (comparative, eventalign)
│   ├── nanocompore (comparative, eventalign)
│   ├── differr (comparative)
│   ├── drummer (comparative)
│   ├── eligos2 (comparative)
│   └── epinano (comparative)
├── Post-processing (format, annotate)
└── QC (quantification, coverage, variants)
```

### Target Architecture

```
NanoRNAMod Workflow (Expanded)
├── Data Preparation
├── Eventalign (f5c) - For: baleen, nanocompore
├── Basecalling Alternatives - For: Dorado, DirectRM, RNANO, TandemMod
├── Modification Detection (17 tools)
│   ├── Comparative (8 tools)
│   │   ├── xpore ✓
│   │   ├── baleen ✓
│   │   ├── nanocompore ✓
│   │   ├── differr ✓
│   │   ├── drummer ✓
│   │   ├── eligos2 ✓
│   │   ├── epinano ✓
│   │   └── PsiNanopore (NEW)
│   └── Single-Sample (9 tools)
│       ├── Dorado (NEW)
│       ├── DirectRM (NEW)
│       ├── TandemMod (NEW)
│       ├── m6ATM (NEW)
│       ├── RNANO (NEW)
│       ├── NanoPsu (NEW)
│       ├── NanoMUD (NEW)
│       └── Penguin (NEW)
├── Post-processing (format, annotate)
└── QC
```

### Key Design Decisions

1. **Dual Data Pathways**: Support both eventalign-dependent (baleen, nanocompore) and raw signal-based tools (Dorado, DirectRM, etc.)
2. **Configuration Flexibility**: Allow users to enable/disable each tool independently
3. **Backward Compatibility**: Maintain existing functionality while adding new tools
4. **Standardized Output**: Ensure all tools produce consistent output format (TSV with standard columns)

---

## Phase-by-Phase Implementation Plan

### Phase 1: Multi-Modification Tools (Weeks 1-3)

**Objective**: Add state-of-the-art multi-modification detection capabilities

**Tools**:
1. **TandemMod** (Nature 2024, 67 citations)
2. **DirectRM** (Nature 2025, 3 citations)

**Priority**: Highest - These tools detect 3-6 modifications each

**Deliverables**:
- [ ] `workflow/rules/modetect_tandemmod.smk`
- [ ] `workflow/rules/modetect_directrm.smk`
- [ ] `workflow/envs/tandemmod.yaml`
- [ ] `workflow/envs/directrm.yaml`
- [ ] `workflow/scripts/tandemmod_*.py` (helper scripts)
- [ ] `workflow/scripts/directrm_*.py` (helper scripts)
- [ ] Update `config/config.yaml` with tool configurations
- [ ] Update `workflow/Snakefile` includes
- [ ] Update `workflow/rules/common.smk` get_final_output()
- [ ] Update `workflow/rules/post_format.smk` formatting rules

**Technical Requirements**:
- TandemMod: Python 3.8+, torch, pandas, numpy, scikit-learn
- DirectRM: Python 3.8+, deep learning framework (PyTorch/TensorFlow), pandas
- Both: Input from BAM files (no eventalign required)
- Both: Output modification predictions per transcript/position

**Implementation Steps**:
1. Research tool installation and usage (GitHub repos, documentation)
2. Create conda environment files with dependencies
3. Develop rule files:
   - Input: BAM files from minimap2 alignment
   - Process: Run tool on each sample (single-sample mode)
   - Output: Modification predictions
4. Create helper scripts for configuration and execution
5. Add formatting rules in post_format.smk
6. Update common.smk to include in final outputs
7. Test with sample data

**Testing**:
- Unit tests for each rule
- Integration test with existing workflow
- Validate output format consistency
- Benchmark resource usage

---

### Phase 2: Enhanced m6A Detection (Weeks 4-5)

**Objective**: Improve m6A detection accuracy with modern deep learning approaches

**Tools**:
3. **m6ATM** (2024, 11 citations) - Deep learning, better than m6Anet
4. **RNANO** (bioRxiv 2025) - NLP-based approach

**Priority**: High - Significant improvement over existing m6A tools

**Deliverables**:
- [ ] `workflow/rules/modetect_m6atm.smk`
- [ ] `workflow/rules/modetect_rnano.smk`
- [ ] `workflow/envs/m6atm.yaml`
- [ ] `workflow/envs/rnano.yaml`
- [ ] `workflow/scripts/m6atm_*.py`
- [ ] `workflow/scripts/rnano_*.py`
- [ ] Config updates
- [ ] Common.smk updates
- [ ] Post-format updates

**Technical Requirements**:
- m6ATM: Python 3.8+, WaveNet, noisy-OR pooling, pandas, numpy
- RNANO: Python 3.9+, transformers (BERT/ALBERT), pandas, huggingface
- Both: Input from BAM files
- Both: Output m6A predictions with confidence scores

**Implementation Steps**:
1. Research m6ATM GitHub repo and installation
2. Research RNANO implementation and models
3. Create conda environments (m6ATM: torch, RNANO: transformers)
4. Develop rule files for single-sample analysis
5. Create helper scripts for tool execution
6. Add formatting and annotation rules
7. Integration testing

**Special Considerations**:
- RNANO may require GPU for optimal performance (optional, CPU fallback)
- Both tools provide stoichiometry estimates - add to output format
- Compare performance with existing m6A tools (xpore, epinano)

---

### Phase 3: Pseudouridine Detection (Weeks 6-8)

**Objective**: Comprehensive Ψ detection with multiple specialized tools

**Tools**:
5. **PsiNanopore** (Nature 2023, 120 citations) - Highest Ψ accuracy
6. **NanoPsu** (Genome Biology 2021, 118 citations) - Standard Ψ tool
7. **NanoMUD** (2024, 8 citations) - Ψ + N1-methylpseudouridine (m1Ψ)
8. **Penguin** (Methods 2022, 57 citations) - Popular Ψ tool

**Priority**: High - 4 specialized Ψ tools provide comprehensive coverage

**Deliverables**:
- [ ] `workflow/rules/modetect_psipore.smk` (PsiNanopore)
- [ ] `workflow/rules/modetect_nanopsu.smk`
- [ ] `workflow/rules/modetect_nanomud.smk`
- [ ] `workflow/rules/modetect_penguin.smk`
- [ ] 4 conda environment files
- [ ] 8-12 helper scripts
- [ ] Config updates
- [ ] Common.smk updates
- [ ] Post-format updates

**Tool-Specific Requirements**:

**PsiNanopore**:
- Type: Comparative (requires control samples)
- Input: Native/Control BAM files
- Dependencies: Python 3.7+, numpy, pandas, scikit-learn
- Output: Ψ predictions with p-values

**NanoPsu**:
- Type: Single-sample
- Input: BAM files
- Dependencies: Python 3.7+, keras/tensorflow, pandas
- Output: Ψ predictions with stoichiometry

**NanoMUD**:
- Type: Single-sample
- Input: BAM files
- Dependencies: Python 3.8+, torch, pandas, numpy
- Output: Ψ + m1Ψ predictions (unique modification)
- Features: AUC 0.998, exceptional accuracy

**Penguin**:
- Type: Single-sample
- Input: BAM files
- Dependencies: Python 3.7+, xgboost/lightgbm, pandas
- Output: Ψ predictions
- Features: Ensemble ML models

**Implementation Steps**:
1. Research each tool's GitHub repo and documentation
2. Create conda environments (PsiNanopore: sklearn, NanoPsu: tensorflow, NanoMUD: torch, Penguin: xgboost)
3. Develop rule files:
   - PsiNanopore: Comparative (native/control)
   - Others: Single-sample
4. Create helper scripts for configuration and execution
5. Add NanoMUD special handling for m1Ψ modification
6. Add formatting rules for each tool's output format
7. Update common.smk for all 4 tools
8. Comprehensive testing

**Special Considerations**:
- NanoMUD is only tool detecting m1Ψ - highlight in documentation
- PsiNanopore requires negative control samples - ensure workflow supports this
- Different ML frameworks across tools (test compatibility)

---

### Phase 4: ONT Native Integration & Additional Tools (Weeks 9-10)

**Objective**: Add ONT's native tool and m6A accuracy improvements

**Tools**:
9. **Dorado** (ONT 2024) - Built-in basecaller with modification detection
10. **EpiNano** (2021, 87 citations) - High-accuracy m6A

**Priority**: Medium - ONT native integration, complementary m6A tool

**Deliverables**:
- [ ] `workflow/rules/modetect_dorado.smk`
- [ ] `workflow/rules/modetect_epinano2.smk` (rename to avoid conflict)
- [ ] `workflow/envs/dorado.yaml`
- [ ] `workflow/envs/epinano2.yaml`
- [ ] Helper scripts
- [ ] Config updates
- [ ] Common.smk updates
- [ ] Post-format updates

**Tool-Specific Requirements**:

**Dorado**:
- Type: Single-sample, built-in basecaller
- Input: FAST5 files (no pre-alignment needed)
- Dependencies: Python 3.8+, ONT Dorado package
- Output: Basecalled sequences with modification probabilities
- Chemistry: RNA004 (new) and RNA002 (legacy)
- Features: Real-time modification detection, no separate installation

**EpiNano**:
- Type: Comparative
- Input: Native/Control BAM files
- Dependencies: Python 3.7+, sklearn, pandas
- Output: m6A predictions (~90% accuracy)

**Implementation Steps**:
1. Dorado: Research installation and usage patterns
   - Option A: Replace guppy with Dorado for basecalling
   - Option B: Use Dorado on pre-aligned data (if supported)
   - Recommended: Option A for full RNA004 support
2. Create conda environment for Dorado
3. Create rule file for Dorado basecalling + modification detection
4. EpiNano: Research GitHub repo
5. Create conda environment
6. Create rule file for comparative analysis
7. Add formatting rules
8. Integration and testing

**Special Considerations**:
- Dorado may replace guppy basecalling - consider workflow restructuring
- Dorado supports RNA004 chemistry natively - document as key feature
- EpiNano provides alternative m6A method - good for comparison

---

## Tool-Specific Implementation Details

### Tool Implementation Template

For each tool, implement:

```yaml
# 1. Rule File Structure (workflow/rules/modetect_{tool}.smk)
rule {tool}_run:
    input:
        # Tool-specific inputs
    output:
        # Standardized output
    params:
        # Configuration parameters
    threads: config["threads"]["{tool}"]
    resources:
        mem_mb = <memory_requirement>
    log:
        "logs/{project}/{tool}/{sample}.log"
    conda:
        "../envs/{tool}.yaml"
    shell/script:
        # Tool execution
```

```yaml
# 2. Conda Environment (workflow/envs/{tool}.yaml)
channels:
  - conda-forge
  - bioconda
dependencies:
  - tool_name = version
  - dependency1
  - dependency2
```

```yaml
# 3. Config Entry (config/config.yaml)
tools:
  {tool}:
    activate: True
params:
  {tool}: "parameter_string"
threads:
  {tool}: 4
```

```python
# 4. Format Rule (workflow/rules/post_format.smk)
rule post_{tool}:
    input:
        "{project}/results/{tool}/{native}_{control}/tool_output"
    output:
        "{project}/results/modifications/{tool}/{native}_{control}/{tool}_results.tsv"
    params:
        tool="{tool}"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"
```

### Detailed Implementation for Each Tool

#### 1. TandemMod

**Repository**: Search GitHub for "TandemMod nanopore"
**Installation**: `pip install tandemmod` or clone from repo
**Input**: BAM files (aligned reads)
**Command**: `tandemmod predict --input bam --output results --model pretrained`
**Output Format**: TSV with columns (transcript, position, modification, probability)
**Special Requirements**: Pre-trained models download

**Implementation**:
```smk
# workflow/rules/modetect_tandemmod.smk
rule tandemmod_predict:
    input:
        bam="{project}/results/bam/{sample}.bam",
        reference=config["reference"]["transcriptome_fasta"]
    output:
        directory="{project}/results/tandemmod/{sample}_predictions")
    threads: config["threads"]["tandemmod"]
    resources:
        mem_mb = 1024 * 100
    params:
        extra="",
        model_type="multi"  # m6A, Ψ, m5C
    conda:
        "../envs/tandemmod.yaml"
    shell:
        "tandemmod predict "
        "--bam {input.bam} "
        "--reference {input.reference} "
        "--model {params.model_type} "
        "--output {output.directory} "
        "--threads {threads} "
        "{params.extra} "
        "2>{log}"
```

#### 2. DirectRM

**Repository**: Search GitHub for "DirectRM nanopore"
**Installation**: Clone and install dependencies
**Input**: BAM files
**Command**: `directrm predict --bam input.bam --output results --mods all`
**Output Format**: TSV with 6 modifications (m6A, Ψ, m5C, A-to-I, m7G, m1A)
**Special Requirements**: Large model files (download on first run)

**Implementation**:
```smk
# workflow/rules/modetect_directrm.smk
rule directrm_predict:
    input:
        bam="{project}/results/bam/{sample}.bam",
        reference=config["reference"]["transcriptome_fasta"]
    output:
        directory="{project}/results/directrm/{sample}_predictions")
    threads: config["threads"]["directrm"]
    resources:
        mem_mb = 1024 * 100
    params:
        modifications="all",  # 6 modifications
        extra=""
    conda:
        "../envs/directrm.yaml"
    shell:
        "directrm predict "
        "--bam {input.bam} "
        "--reference {input.reference} "
        "--mods {params.modifications} "
        "--output {output.directory} "
        "--threads {threads} "
        "{params.extra} "
        "2>{log}"
```

#### 3. m6ATM

**Repository**: Search GitHub for "m6ATM nanopore"
**Installation**: `pip install m6atm` or clone
**Input**: BAM files
**Command**: `m6atm predict --bam input.bam --output results`
**Output Format**: TSV with m6A predictions, stoichiometry
**Special Requirements**: WaveNet encoder, noisy-OR pooling

**Implementation**:
```smk
# workflow/rules/modetect_m6atm.smk
rule m6atm_predict:
    input:
        bam="{project}/results/bam/{sample}.bam"
    output:
        directory="{project}/results/m6atm/{sample}_predictions")
    threads: config["threads"]["m6atm"]
    resources:
        mem_mb = 1024 * 80
    params:
        stoichiometry=True
    conda:
        "../envs/m6atm.yaml"
    shell:
        "m6atm predict "
        "--bam {input.bam} "
        "--output {output.directory} "
        "--stoichiometry {params.stoichiometry} "
        "--threads {threads} "
        "2>{log}"
```

#### 4. RNANO

**Repository**: Search bioRxiv paper for GitHub link
**Installation**: Clone and install
**Input**: BAM files
**Command**: `rnano predict --bam input.bam --output results`
**Output Format**: TSV with multi-mod predictions
**Special Requirements**: Transformers (BERT/ALBERT), GPU optional

**Implementation**:
```smk
# workflow/rules/modetect_rnano.smk
rule rnano_predict:
    input:
        bam="{project}/results/bam/{sample}.bam",
        reference=config["reference"]["transcriptome_fasta"]
    output:
        directory="{project}/results/rnano/{sample}_predictions")
    threads: config["threads"]["rnano"]
    resources:
        mem_mb = 1024 * 100
    params:
        gpu=False,  # Set to True if GPU available
        model="pretrained"
    conda:
        "../envs/rnano.yaml"
    shell:
        "rnano predict "
        "--bam {input.bam} "
        "--reference {input.reference} "
        "--model {params.model} "
        "--gpu {params.gpu} "
        "--output {output.directory} "
        "--threads {threads} "
        "2>{log}"
```

#### 5. PsiNanopore

**Repository**: Search GitHub for "PsiNanopore"
**Installation**: Clone and install
**Input**: Native/Control BAM files
**Command**: `psipore detect --native native.bam --control control.bam --output results`
**Output Format**: TSV with Ψ predictions, p-values
**Special Requirements**: Negative control samples

**Implementation**:
```smk
# workflow/rules/modetect_psipore.smk
rule psipore_detect:
    input:
        native_bam="{project}/results/bam/{native}.bam",
        control_bam="{project}/results/bam/{control}.bam"
    output:
        directory="{project}/results/psipore/{native}_{control}_predictions")
    threads: config["threads"]["psipore"]
    resources:
        mem_mb = 1024 * 50
    params:
        min_coverage=10
    conda:
        "../envs/psipore.yaml"
    shell:
        "psipore detect "
        "--native {input.native_bam} "
        "--control {input.control_bam} "
        "--output {output.directory} "
        "--min_coverage {params.min_coverage} "
        "--threads {threads} "
        "2>{log}"
```

#### 6. NanoPsu

**Repository**: Search GitHub for "NanoPsu"
**Installation**: Clone and install
**Input**: BAM files
**Command**: `nanopsu predict --bam input.bam --output results`
**Output Format**: TSV with Ψ predictions, stoichiometry
**Special Requirements**: Keras/TensorFlow

**Implementation**:
```smk
# workflow/rules/modetect_nanopsu.smk
rule nanopsu_predict:
    input:
        bam="{project}/results/bam/{sample}.bam"
    output:
        directory="{project}/results/nanopsu/{sample}_predictions")
    threads: config["threads"]["nanopsu"]
    resources:
        mem_mb = 1024 * 60
    params:
        stoichiometry=True
    conda:
        "../envs/nanopsu.yaml"
    shell:
        "nanopsu predict "
        "--bam {input.bam} "
        "--output {output.directory} "
        "--stoichiometry {params.stoichiometry} "
        "--threads {threads} "
        "2>{log}"
```

#### 7. NanoMUD

**Repository**: Search GitHub for "NanoMUD"
**Installation**: Clone and install
**Input**: BAM files
**Command**: `nanomud predict --bam input.bam --output results`
**Output Format**: TSV with Ψ and m1Ψ predictions
**Special Requirements**: Only tool detecting m1Ψ

**Implementation**:
```smk
# workflow/rules/modetect_nanomud.smk
rule nanomud_predict:
    input:
        bam="{project}/results/bam/{sample}.bam"
    output:
        directory="{project}/results/nanomud/{sample}_predictions")
    threads: config["threads"]["nanomud"]
    resources:
        mem_mb = 1024 * 80
    params:
        mods=["Psi", "m1Psi"]  # Both modifications
    conda:
        "../envs/nanomud.yaml"
    shell:
        "nanomud predict "
        "--bam {input.bam} "
        "--output {output.directory} "
        "--mods {params.mods} "
        "--threads {threads} "
        "2>{log}"
```

#### 8. Penguin

**Repository**: Search GitHub for "Penguin pseudouridine nanopore"
**Installation**: `pip install penguin-nanopore` or clone
**Input**: BAM files
**Command**: `penguin predict --bam input.bam --output results`
**Output Format**: TSV with Ψ predictions
**Special Requirements**: XGBoost/LightGBM

**Implementation**:
```smk
# workflow/rules/modetect_penguin.smk
rule penguin_predict:
    input:
        bam="{project}/results/bam/{sample}.bam"
    output:
        directory="{project}/results/penguin/{sample}_predictions")
    threads: config["threads"]["penguin"]
    resources:
        mem_mb = 1024 * 50
    params:
        model="ensemble"
    conda:
        "../envs/penguin.yaml"
    shell:
        "penguin predict "
        "--bam {input.bam} "
        "--output {output.directory} "
        "--model {params.model} "
        "--threads {threads} "
        "2>{log}"
```

#### 9. Dorado

**Repository**: ONT official or GitHub: `https://github.com/nanoporetech/dorado`
**Installation**: `conda install -c bioconda dorado` or download binary
**Input**: FAST5 files (raw signal)
**Command**: `dorado basecaller fast5s/ model_ref/ output/ --modified-bases model`
**Output Format**: BAM file with modification tags
**Special Requirements**: RNA004 chemistry support, model files download

**Implementation Options**:

**Option A: Replace Guppy (Recommended for RNA004)**:
```smk
# workflow/rules/prep_basecall_dorado.smk (NEW FILE)
rule dorado_basecall:
    input:
        fast5="{project}/data/{sample}/fast5/",
        reference=config["reference"]["transcriptome_fasta"]
    output:
        bam="{project}/results/bam/{sample}.bam",
        bai="{project}/results/bam/{sample}.bam.bai"
    threads: config["threads"]["dorado"]
    resources:
        mem_mb = 1024 * 50
    params:
        model="rna004_130bps_hac@v5.0.0",  # RNA004 model
        modified_bases=True
    conda:
        "../envs/dorado.yaml"
    shell:
        "dorado basecaller "
        "{input.fast5} "
        "{params.model} "
        "--modified-bases {params.modified_bases} "
        "--reference {input.reference} "
        "--output-dir {output.bam} "
        "--threads {threads} "
        "&& samtools index {output.bam} "
        "2>{log}"
```

**Option B: Use with pre-aligned BAM (if supported)**:
```smk
# workflow/rules/modetect_dorado.smk
rule dorado_modcall:
    input:
        bam="{project}/results/bam/{sample}.bam",
        reference=config["reference"]["transcriptome_fasta"]
    output:
        directory="{project}/results/dorado/{sample}_modifications")
    threads: config["threads"]["dorado"]
    resources:
        mem_mb = 1024 * 60
    params:
        model="rna004_130bps_hac@v5.0.0"
    conda:
        "../envs/dorado.yaml"
    shell:
        "dorado modifier "
        "--bam {input.bam} "
        "--model {params.model} "
        "--output {output.directory} "
        "--threads {threads} "
        "2>{log}"
```

**Recommendation**: Implement Option A with a config flag to choose between guppy and dorado basecalling.

#### 10. EpiNano (rename to epinano2)

**Repository**: Search GitHub for "EpiNano"
**Installation**: Clone and install
**Input**: Native/Control BAM files
**Command**: `epinano detect --native native.bam --control control.bam --output results`
**Output Format**: TSV with m6A predictions
**Special Requirements**: High accuracy (~90%)

**Implementation**:
```smk
# workflow/rules/modetect_epinano2.smk
rule epinano2_detect:
    input:
        native_bam="{project}/results/bam/{native}.bam",
        control_bam="{project}/results/bam/{control}.bam"
    output:
        directory="{project}/results/epinano2/{native}_{control}_predictions")
    threads: config["threads"]["epinano2"]
    resources:
        mem_mb = 1024 * 50
    params:
        min_reads=100
    conda:
        "../envs/epinano2.yaml"
    shell:
        "epinano detect "
        "--native {input.native_bam} "
        "--control {input.control_bam} "
        "--output {output.directory} "
        "--min_reads {params.min_reads} "
        "--threads {threads} "
        "2>{log}"
```

---

## Technical Specifications

### Standard Output Format

All tools must produce output in standardized TSV format:

```tsv
transcript	position	modification	probability	score	stoichiometry	other_columns...
chr1	123	m6A	0.95	0.45	...
chr1	456	Ψ	0.87	0.32	...
```

**Minimum Required Columns**:
- `transcript`: Transcript ID
- `position`: Modification position (0-based)
- `modification`: Modification type (m6A, Ψ, m5C, etc.)
- `probability` or `score`: Confidence value

**Optional Columns**:
- `stoichiometry`: Modification fraction
- `p_value`: Statistical significance (comparative tools)
- `strand`: + or -
- `context`: Sequence context

### Memory and Thread Requirements

| Tool | Memory (GB) | Threads | Recommended | Notes |
|-------|--------------|---------|--------------|--------|
| TandemMod | 100 | 4-8 | High | Deep learning |
| DirectRM | 100 | 4-8 | High | 6 modifications |
| m6ATM | 80 | 4-8 | High | WaveNet |
| RNANO | 100 | 4-8 | High | Transformers |
| PsiNanopore | 50 | 4-8 | Medium | Comparative |
| NanoPsu | 60 | 4-8 | Medium | TensorFlow |
| NanoMUD | 80 | 4-8 | Medium | 2 modifications |
| Penguin | 50 | 4-8 | Medium | XGBoost |
| Dorado | 50-100 | 4-16 | High | Basecalling |
| EpiNano2 | 50 | 4-8 | Medium | Comparative |

### Dependency Conflicts

**Potential Conflicts**:
- TensorFlow vs PyTorch: NanoPsu (TF) vs others (PyTorch)
- Multiple ML frameworks in same environment
- GPU requirements (RNANO, DirectRM, m6ATM)

**Resolution Strategy**:
- Separate conda environments for each tool (current approach)
- CPU-only versions as fallback
- Document GPU requirements clearly

---

## Testing and Validation Strategy

### Unit Testing

For each new tool:

1. **Installation Test**:
   ```bash
   conda create -n test_{tool} --file workflow/envs/{tool}.yaml
   conda activate test_{tool}
   python -c "import {tool}_package"
   ```

2. **Execution Test**:
   ```bash
   snakemake --dry-run -R {tool}_predict
   snakemake -R {tool}_predict --cores 1
   ```

3. **Output Validation**:
   - Check file exists
   - Validate TSV format
   - Verify required columns present
   - Check no NaN/Inf values

### Integration Testing

1. **Full Workflow Test**:
   ```bash
   # Enable single tool
   cd .test
   snakemake --use-conda --cores 4 -k
   ```

2. **Multi-Tool Test**:
   ```bash
   # Enable multiple tools
   # Run full workflow
   # Check all outputs present
   ```

3. **Comparative Test**:
   - Run same dataset with multiple tools
   - Compare results
   - Check consistency

### Benchmark Testing

1. **Resource Usage**:
   - Monitor CPU, memory, disk I/O
   - Compare to current tools
   - Optimize if necessary

2. **Runtime**:
   - Measure execution time
   - Compare to literature benchmarks
   - Document expected runtime

3. **Accuracy**:
   - Test on simulated data with known modifications
   - Compare to ground truth
   - Document accuracy metrics

---

## Dependency Management

### Conda Environment Files

Create 10 new environment files:

```
workflow/envs/
├── tandemmod.yaml (NEW)
├── directrm.yaml (NEW)
├── m6atm.yaml (NEW)
├── rnano.yaml (NEW)
├── psipore.yaml (NEW)
├── nanopsu.yaml (NEW)
├── nanomud.yaml (NEW)
├── penguin.yaml (NEW)
├── dorado.yaml (NEW)
└── epinano2.yaml (NEW)
```

**Template**:
```yaml
name: {tool}
channels:
  - conda-forge
  - bioconda
  - nodefaults
dependencies:
  - python=3.9
  - {tool_name}
  - pandas
  - numpy
  - {other_dependencies}
  - pip:
      - {pip_only_packages}
```

### Script Dependencies

**Helper Scripts to Create**:
- `tandemmod_config.py` - Generate TandemMod config files
- `directrm_config.py` - Generate DirectRM config files
- `m6atm_predict.py` - Wrapper for m6ATM execution
- `rnano_predict.py` - Wrapper for RNANO execution
- `psipore_config.py` - Generate PsiNanopore config
- `nanopsu_predict.py` - Wrapper for NanoPsu
- `nanomud_predict.py` - Wrapper for NanoMUD
- `penguin_predict.py` - Wrapper for Penguin
- `dorado_download_models.py` - Download Dorado models
- `epinano2_config.py` - Generate EpiNano2 config

### External Tool Installation

**Download and Install Scripts**:
- Create `workflow/scripts/download_tools.py`
- Automatically download tools from GitHub
- Install to appropriate locations
- Handle version pinning

---

## Configuration Updates

### config.yaml Updates

Add to existing `config.yaml`:

```yaml
# New tool configurations
params:
  # Existing params...
  tandemmod: "--model multi --threshold 0.5"
  directrm: "--mods all --min-prob 0.7"
  m6atm: "--stoichiometry --threshold 0.6"
  rnano: "--model pretrained --gpu False"
  psipore: "--min-coverage 10 --p-value 0.05"
  nanopsu: "--stoichiometry --min-prob 0.5"
  nanomud: "--mods Psi,m1Psi --threshold 0.8"
  penguin: "--model ensemble --threshold 0.6"
  dorado: "--model rna004_130bps_hac@v5.0.0 --modified-bases"
  epinano2: "--min-reads 100 --p-value 0.05"

threads:
  # Existing threads...
  tandemmod: 4
  directrm: 4
  m6atm: 4
  rnano: 4
  psipore: 4
  nanopsu: 4
  nanomud: 4
  penguin: 4
  dorado: 8  # Higher for basecalling
  epinano2: 4

tools:
  # Existing tools...
  tandemmod:
    activate: False  # Default off, user enables
  directrm:
    activate: False
  m6atm:
    activate: False
  rnano:
    activate: False
  psipore:
    activate: False
  nanopsu:
    activate: False
  nanomud:
    activate: False
  penguin:
    activate: False
  dorado:
    activate: False
    basecaller: False  # Option to use Dorado for basecalling
  epinano2:
    activate: False
```

### workflow/Snakefile Updates

Add includes:

```python
# Add after existing includes
include: "rules/modetect_tandemmod.smk"
include: "rules/modetect_directrm.smk"
include: "rules/modetect_m6atm.smk"
include: "rules/modetect_rnano.smk"
include: "rules/modetect_psipore.smk"
include: "rules/modetect_nanopsu.smk"
include: "rules/modetect_nanomud.smk"
include: "rules/modetect_penguin.smk"
include: "rules/modetect_dorado.smk"
include: "rules/modetect_epinano2.smk"
```

**Optional: Add new basecalling rule** (if using Dorado):

```python
include: "rules/prep_basecall_dorado.smk"  # NEW FILE
```

### workflow/rules/common.smk Updates

Add to `get_final_output()` function:

```python
def get_final_output():
    # Existing code...
    tools = [tool for tool in config["tools"] if config["tools"][tool]["activate"]]

    # Add new tools
    if "tandemmod" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/tandemmod/{sample}/tandemmod_results.tsv",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/tandemmod/{sample}/tandemmod_annotated_results.tsv",
                sample=list(samples.index),
                RESULT_ROOT=RESULT_ROOT,
            )

    if "directrm" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/directrm/{sample}/directrm_results.tsv",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/directrm/{sample}/directrm_annotated_results.tsv",
                sample=list(samples.index),
                RESULT_ROOT=RESULT_ROOT,
            )

    # Repeat for all 10 tools...

    return final_output
```

### workflow/rules/post_format.smk Updates

Add formatting rules for each tool:

```smk
# Example for single-sample tools
rule post_tandemmod:
    input:
        "{project}/results/tandemmod/{sample}_predictions/predictions.tsv",
        "{project}/results/tandemmod/{sample}_predictions"
    output:
        "{project}/results/modifications/tandemmod/{sample}/tandemmod_results.tsv",
    params:
        tool="tandemmod",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule tandemmod_annotate:
    input:
        transcripts="{project}/results/modifications/tandemmod/{sample}/tandemmod_results.tsv",
    output:
        result="{project}/results/modifications/tandemmod/{sample}/tandemmod_annotated_results.tsv",
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} "
        "2>logs/{project}/tandemmod_annotate/{sample}.err "
        "1>logs/{project}/tandemmod_annotate/{sample}.log"

# Repeat for all 10 tools...
```

**Update format.py**:

```python
# Add new tool formatting functions
def format_tandemmod(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    df = df.sort_values(['transcript', 'position'])
    df.to_csv(output_file, sep='\t', index=False)

def format_directrm(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    df = df.sort_values(['transcript', 'position'])
    df.to_csv(output_file, sep='\t', index=False)

# Add for all 10 tools...

if snakemake.params.tool == 'tandemmod':
    format_tandemmod(snakemake.input[0], snakemake.output[0])
# Add for all 10 tools...
```

---

## Documentation Requirements

### AGENTS.md Updates

Update `AGENTS.md` with new tools:

```markdown
### Nanopore RNA Modification Tools

**Currently Supported (17 tools)**:

**Comparative Tools**:
- xpore, baleen, nanocompore, differr, drummer, eligos2, epinano, epinano2, psipore

**Single-Sample Tools**:
- tandemmod, directrm, m6atm, rnano, nanopsu, nanomud, penguin

**Basecalling with Modification Detection**:
- dorado (RNA004 chemistry support)

**Tool Selection**:
- Enable tools in config.yaml under `tools` section
- Each tool can be enabled/disabled independently
- All tools produce standardized TSV output format
```

### README.md Updates

Add section:

```markdown
## Supported Tools

The workflow supports 17 Nanopore-based RNA modification detection tools:

### Multi-Modification Tools
- **TandemMod**: Detects m6A, Ψ, m5C using transfer learning
- **DirectRM**: Detects 6 modifications (m6A, Ψ, m5C, A-to-I, m7G, m1A)

### Specialized m6A Tools
- **m6ATM**: Deep learning with WaveNet encoder (better than m6Anet)
- **RNANO**: NLP-based approach using transformers
- **EpiNano2**: High-accuracy comparative analysis

### Pseudouridine (Ψ) Tools
- **PsiNanopore**: Highest accuracy Ψ detection (comparative)
- **NanoPsu**: Standard Ψ detection with stoichiometry
- **NanoMUD**: Unique - detects Ψ and N1-methylpseudouridine (m1Ψ)
- **Penguin**: Ensemble ML models for Ψ detection

### Basecalling Integration
- **Dorado**: ONT native basecaller with modification detection (RNA004 support)

### Existing Tools
- xpore, baleen, nanocompore, differr, drummer, eligos2, epinano

## Tool Activation

To enable a tool, set `activate: True` in config/config.yaml:

```yaml
tools:
  tandemmod:
    activate: True
  dorado:
    activate: True
    basecaller: False  # Use Dorado for basecalling
```

## RNA004 Chemistry Support

Dorado supports the new RNA004 chemistry (released 2024):
- Requires RNA-specific flowcell
- Native modification detection
- Better accuracy for many modification types

To enable:
```yaml
tools:
  dorado:
    activate: True
    basecaller: True  # Use Dorado instead of Guppy
params:
  dorado: "--model rna004_130bps_hac@v5.0.0"
```
```

### Tool-Specific Documentation

Create `workflow/docs/TOOLS.md`:

```markdown
# Tool-Specific Documentation

## TandemMod
### Publication
Wu et al., Nature 2024

### Installation
```bash
pip install tandemmod
```

### Usage in NanoRNAMod
```yaml
tools:
  tandemmod:
    activate: True
params:
  tandemmod: "--model multi --threshold 0.5"
```

### Output Format
- transcript: Transcript ID
- position: Modification position (0-based)
- modification: m6A, Ψ, m5C
- probability: Confidence score (0-1)
- stoichiometry: Modification fraction (0-1)

### Parameters
- `--model`: multi (all 3 modifications) or individual
- `--threshold`: Minimum probability (default 0.5)
- `--threads`: Number of threads

[Repeat for all 10 tools...]
```

---

## Timeline and Milestones

### Phase 1: Multi-Modification Tools (Weeks 1-3)

**Week 1**:
- [ ] Research TandemMod and DirectRM (installation, usage)
- [ ] Create conda environment files
- [ ] Create helper scripts (config, predict)
- [ ] Create rule files for TandemMod

**Week 2**:
- [ ] Create rule files for DirectRM
- [ ] Implement formatting rules
- [ ] Update common.smk get_final_output()
- [ ] Unit testing

**Week 3**:
- [ ] Integration testing
- [ ] Benchmark testing
- [ ] Documentation updates
- [ ] Phase 1 review and approval

**Milestone**: 2 new multi-modification tools operational

### Phase 2: Enhanced m6A Detection (Weeks 4-5)

**Week 4**:
- [ ] Research m6ATM and RNANO
- [ ] Create conda environment files
- [ ] Create helper scripts
- [ ] Create rule files

**Week 5**:
- [ ] Implement formatting rules
- [ ] Update common.smk
- [ ] Testing and validation
- [ ] Documentation updates
- [ ] Phase 2 review and approval

**Milestone**: 2 enhanced m6A tools operational

### Phase 3: Pseudouridine Detection (Weeks 6-8)

**Week 6**:
- [ ] Research all 4 Ψ tools
- [ ] Create conda environment files
- [ ] Create helper scripts for PsiNanopore, NanoPsu

**Week 7**:
- [ ] Create helper scripts for NanoMUD, Penguin
- [ ] Create rule files for all 4 tools
- [ ] Implement formatting rules

**Week 8**:
- [ ] Update common.smk
- [ ] Comprehensive testing
- [ ] Benchmarking
- [ ] Documentation
- [ ] Phase 3 review and approval

**Milestone**: 4 specialized Ψ tools operational

### Phase 4: ONT Integration & Additional Tools (Weeks 9-10)

**Week 9**:
- [ ] Research Dorado (installation, RNA004 support)
- [ ] Research EpiNano2
- [ ] Create conda environment files
- [ ] Create Dorado basecalling rule (optional)

**Week 10**:
- [ ] Create Dorado modcalling rule
- [ ] Create EpiNano2 rule
- [ ] Implement formatting rules
- [ ] Update common.smk
- [ ] Full workflow testing
- [ ] Final documentation
- [ ] Phase 4 review and approval

**Milestone**: Dorado integration and EpiNano2 operational

### Final Review (Week 11-12)

**Week 11**:
- [ ] Full integration testing (all 10 tools)
- [ ] Performance benchmarking
- [ ] Memory and resource optimization
- [ ] Bug fixes

**Week 12**:
- [ ] Final documentation review
- [ ] User guide updates
- [ ] AGENTS.md updates
- [ ] README.md updates
- [ ] Prepare release

**Milestone**: NanoRNAMod 2.0 release with 17 tools

---

## Risk Assessment and Mitigation

### High-Risk Items

| Risk | Probability | Impact | Mitigation Strategy |
|-------|-------------|--------|-------------------|
| Tool repositories not available | Medium | High | Search alternatives, contact authors, create fallback |
| Conflicting dependencies (TF vs PyTorch) | Medium | High | Separate environments, CPU-only versions |
| Insufficient memory for deep learning | Medium | Medium | Optimize batch processing, add memory flags |
| GPU requirements not met | Low | Medium | CPU-only versions, document requirements |
| Model download failures | Low | Medium | Mirror models, retry logic, cache |

### Medium-Risk Items

| Risk | Probability | Impact | Mitigation Strategy |
|-------|-------------|--------|-------------------|
| Output format incompatibility | Medium | Medium | Create adapter scripts, standardize early |
| Tool updates breaking compatibility | Low | Medium | Version pinning, document tested versions |
| RNA004 chemistry data unavailable | Low | High | Support both RNA002 and RNA004, document requirements |
| Comparative tool control samples missing | Low | High | Validate input requirements, clear error messages |

### Low-Risk Items

| Risk | Probability | Impact | Mitigation Strategy |
|-------|-------------|--------|-------------------|
| Documentation gaps | Low | Low | Thorough review, user feedback |
| Performance regression | Low | Low | Benchmark before/after, optimization |
| Naming conflicts (epinano vs epinano2) | Low | Low | Clear naming, documentation |

### Contingency Plans

**If Tool Installation Fails**:
1. Research alternative tools
2. Skip problematic tool in phase
3. Document as "future work"
4. Continue with remaining tools

**If Integration Testing Fails**:
1. Isolate to specific tool
2. Debug rule file
3. Check conda environment
4. Consult tool documentation

**If Performance is Unacceptable**:
1. Optimize parameters
2. Reduce data scope for testing
3. Add resource limits
4. Document limitations

---

## Success Criteria

### Functional Requirements

- [ ] All 10 tools successfully added to workflow
- [ ] Each tool can be enabled/disabled independently
- [ ] All tools produce standardized output format
- [ ] Output annotation works for all tools
- [ ] Integration with existing rules (prep_data, QC) works

### Non-Functional Requirements

- [ ] Memory usage < specified limits
- [ ] Runtime comparable to similar tools
- [ ] No dependency conflicts
- [ ] Documentation is complete and clear
- [ ] AGENTS.md updated with new tools

### Quality Requirements

- [ ] Unit tests pass for all 10 tools
- [ ] Integration tests pass
- [ ] Benchmark results documented
- [ ] Output validated on test data
- [ ] Code follows project conventions

### User Experience Requirements

- [ ] Config file is intuitive
- [ ] Error messages are clear
- [ ] Documentation explains tool activation
- [ ] Tool comparison is easy
- [ ] Results are easy to interpret

---

## Appendices

### Appendix A: Tool Research Checklist

For each tool, document:
- [ ] GitHub repository URL
- [ ] Installation instructions
- [ ] Required Python version
- [ ] Required dependencies
- [ ] Input file requirements
- [ ] Output file format
- [ ] Command-line interface
- [ ] Documentation links
- [ ] Citation information
- [ ] License information

### Appendix B: Testing Checklist

For each tool, test:
- [ ] Installation works
- [ ] Conda environment creates successfully
- [ ] Rule executes without errors
- [ ] Output file exists
- [ ] Output format is correct
- [ ] Required columns present
- [ ] No NaN/Inf values
- [ ] Integration with post_format works
- [ ] Annotation works if GTF provided
- [ ] Resource usage acceptable

### Appendix C: Code Review Checklist

For each tool's code, review:
- [ ] Follows project naming conventions
- [ ] Uses KEEP_OR_NOT correctly
- [ ] Has proper log files
- [ ] Has benchmark files
- [ ] Resources defined correctly
- [ ] Threads from config
- [ ] Conda environment specified
- [ ] Error handling present
- [ ] Comments are minimal (project convention)
- [ ] No type hints (project convention)

### Appendix D: Release Checklist

Before release:
- [ ] All features implemented
- [ ] All tests pass
- [ ] Documentation complete
- [ ] AGENTS.md updated
- [ ] README.md updated
- [ ] Changelog created
- [ ] Version updated
- [ ] Git commit created
- [ ] Tag created
- [ ] Release notes prepared

---

## Conclusion

This implementation plan provides a **systematic, phased approach** to adding 10 high-priority Nanopore RNA modification tools to the NanoRNAMod workflow. The plan addresses:

✓ **Architecture** - Dual data pathways (eventalign and raw signal)
✓ **Priority** - Multi-modification tools first, then specialized tools
✓ **Technical details** - Specific implementations for each tool
✓ **Testing** - Unit, integration, and benchmark testing
✓ **Documentation** - Comprehensive updates to all docs
✓ **Timeline** - 10-12 weeks with clear milestones
✓ **Risk management** - Identified risks with mitigation strategies

Following this plan will expand NanoRNAMod from 7 to **17 tools**, supporting **8+ modification types** and both **RNA002 and RNA004 chemistries**, significantly enhancing the workflow's capabilities and value to the research community.
