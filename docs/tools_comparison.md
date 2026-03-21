# Nanopore DRS RNA Modification Detection Tools Comparison

This document provides a comprehensive comparison of RNA modification detection tools integrated into NanoRNAMod workflow.

## Overview

Oxford Nanopore direct RNA sequencing (DRS) enables detection of RNA modifications by analyzing changes in ionic current signals and basecalling errors. Tools in this workflow can be categorized by their detection approach and input requirements.

**Key References:**
- [Systematic comparison of tools for m6A mapping (Nature Communications 2023)](https://www.nature.com/articles/s41467-023-37596-5)
- [Comparative evaluation of computational models (Briefings in Bioinformatics 2024)](https://academic.oup.com/bib/article/26/4/bbaf404/8233723)
- [Nanopore DRS for RNA modification analysis (2025)](https://link.springer.com/article/10.1007/s44307-025-00093-5)

---

## Tool Categories

### Category 1: Signal-Based Comparative Tools (Require Native + Control)

These tools compare signal features between modified (Native) and unmodified (Control/IVT) samples.

| Tool | Input | Signal Source | Comparison Type |
|------|-------|---------------|-----------------|
| xPore | eventalign | f5c/nanopolish | GMM-based |
| Nanocompore | eventalign (collapsed) | f5c/nanopolish | Statistical test |
| Baleen | eventalign (indexed) | f5c/nanopolish | HMM + DTW |
| pyBaleen | BAM + BLOW5 + FASTQ | raw signal | CUDA-DTW + HMM |

### Category 2: Error-Based Comparative Tools (Require Native + Control)

These tools analyze basecalling errors/mismatches between samples.

| Tool | Input | Detection Method |
|------|-------|------------------|
| DiffErr | BAM | Basecalling error rate differences |
| DRUMMER | BAM | Basecalling error differences |
| ELIGOS2 | BAM | Error pattern comparison |
| EpiNano | BAM (per-site CSV) | SVM on error features |
| psipore | BAM | Error-based (Psi-specific) |

### Category 3: Single-Sample ML/DL Tools (No Control Required)

These tools use pre-trained models to detect modifications from single samples.

| Tool | Input | Model Type | Modifications |
|------|-------|------------|---------------|
| TandemMod | BAM | Deep Learning (Transfer) | Multi-modification |
| DirectRM | BAM | ML Model | Multiple |
| m6ATM | BAM | Transformer | m6A |
| Rnano | BAM | Pretrained Model | m6A |
| NanoPSU | BAM | ML Model | Pseudouridine (Ψ) |
| NanoMUD | BAM | ML Model | Ψ, m1Ψ |
| Penguin | BAM | ML Model | Pseudouridine (Ψ) |

---

## Detailed Tool Information

### 1. xPore

**Publication:** [Genome Biology 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02593-6)

**Input Requirements:**
- f5c eventalign output (simple format)
- Transcriptome FASTA reference
- Both Native and Control samples

**Output Format:**
```
transcript    position    kmer    statistic    p_value    adjusted_p_value    direction
```

**Supported Modifications:** m6A (primary), other modifications detectable

**Pros:**
- Well-documented and widely used
- GMM-based statistical approach
- Provides per-kmer resolution
- Includes postprocessing for consensus calling

**Cons:**
- Requires eventalign preprocessing (computationally expensive)
- Limited to comparative mode
- Memory-intensive for large datasets

**Benchmarking Considerations:**
- Accuracy depends on eventalign quality
- Sensitive to coverage depth
- kmer-based approach may miss modifications in low-complexity regions

---

### 2. Nanocompore

**Publication:** [Nature Communications 2021](https://www.nature.com/articles/s41467-021-27393-3)

**Input Requirements:**
- f5c eventalign output (full format)
- Transcriptome FASTA reference
- Both Native and Control samples

**Output Format:**
```
transcript    position    ref_kmer    coverage    p_value    p_value_glm    p_value_ks    pass_glm    pass_ks
```

**Supported Modifications:** m6A, can detect other modifications

**Pros:**
- Multiple statistical tests (GLM, KS test)
- Comprehensive QC outputs
- Well-validated on synthetic data

**Cons:**
- Eventalign collapse step is memory-intensive
- Requires significant preprocessing time
- Large intermediate files

**Benchmarking Considerations:**
- Consider both GLM and KS test results
- Coverage thresholds impact sensitivity
- Position-level vs kmer-level resolution

---

### 3. Baleen

**Publication:** [Bioinformatics 2023](https://pubmed.ncbi.nlm.nih.gov/37492785/)

**Input Requirements:**
- f5c eventalign output (bz2 compressed)
- BED file for target regions
- Both Native and Control samples

**Output Format:**
```
transcript    position    pvalue    padj    effect_size    mod_type
```

**Supported Modifications:** m6A and other modifications

**Pros:**
- HMM-based approach with DTW alignment
- Includes postcall refinement step
- Good for transcriptome-wide analysis

**Cons:**
- High memory requirements (450-650GB)
- Multi-step workflow (dataprep, modcall, postcall)
- Complex output structure

**Benchmarking Considerations:**
- Memory usage scales with transcriptome size
- postcall step improves precision
- Target regions BED file affects results

---

### 4. pyBaleen

**GitHub:** [py-baleen](https://github.com/loganylchen/py-baleen)

**Input Requirements:**
- BAM (transcriptome-aligned)
- BLOW5 (raw signal)
- FASTQ
- Both Native and Control samples

**Output Format:**
```
transcript    position    pvalue    padj    effect_size
```

**Supported Modifications:** m6A, general modifications

**Pros:**
- CUDA acceleration for DTW
- Direct raw signal processing
- HMM-based modification calling
- No eventalign preprocessing needed

**Cons:**
- Requires GPU for optimal performance
- High memory requirements
- Experimental/newer tool

**Benchmarking Considerations:**
- CUDA vs CPU performance comparison
- Signal processing quality affects results
- min_depth and min_mapq parameters sensitivity

---

### 5. DiffErr

**Publication:** [Nature Methods 2019](https://www.nature.com/articles/nature2019)

**Input Requirements:**
- Filtered BAM files
- BAM index files
- Transcriptome FASTA
- Both Native and Control samples

**Output Format:** BED format with differential error statistics

**Supported Modifications:** m6A (primary)

**Pros:**
- Simple input requirements (BAM only)
- Fast execution
- Direct error rate comparison

**Cons:**
- Limited to error-based detection
- May miss modifications with subtle error signatures
- Requires matched control

**Benchmarking Considerations:**
- Error rate thresholds affect sensitivity
- Coverage requirements
- Alignment quality critical

---

### 6. DRUMMER

**Publication:** [Bioinformatics 2022](https://pubmed.ncbi.nlm.nih.gov/35426900/)

**GitHub:** [DepledgeLab/DRUMMER](https://github.com/DepledgeLab/DRUMMER)

**Input Requirements:**
- Filtered BAM files
- BAM index files
- Transcriptome FASTA
- Both Native and Control samples

**Output Format:** Summary statistics with modification calls

**Supported Modifications:** m6A, other modifications

**Pros:**
- Nucleotide-level resolution
- Isoform-specific detection
- Fast comparative analysis

**Cons:**
- Requires matched control
- Limited to comparative mode
- Region preparation step needed

**Benchmarking Considerations:**
- Region overlap affects results
- Coverage per isoform
- Statistical threshold selection

---

### 7. ELIGOS2

**Publication:** [Zenodo](https://zenodo.org/records/14199538)

**Input Requirements:**
- Filtered BAM files
- BAM index files
- Transcriptome FASTA
- Merged regions BED file
- Both Native and Control samples

**Output Format:**
```
chr    start    end    ref_base    coverage    error_rate    p_value    q_value
```

**Supported Modifications:** Multiple including m6A, Ψ

**Pros:**
- Can detect various modifications
- Base-level resolution
- Pair-wise comparison mode

**Cons:**
- Requires region preprocessing
- Multiple output files
- Complex parameter tuning

**Benchmarking Considerations:**
- baseExt parameter affects sensitivity
- Region merging strategy impacts results
- Error rate thresholds

---

### 8. EpiNano

**Publication:** [Nature Methods 2019](https://www.nature.com/articles/s41592-019-0617-3)

**GitHub:** [novoalab/EpiNano](https://github.com/novoalab/EpiNano)

**Input Requirements:**
- Per-site error statistics CSV (preprocessed from BAM)
- Reference FASTA with dict file
- Both Native and Control samples

**Output Format:**
```
transcript    position    delta_sum_err    prediction    score
```

**Supported Modifications:** m6A

**Pros:**
- SVM-based classification
- Uses multiple error features
- Well-established method

**Cons:**
- Requires specific preprocessing
- R-based implementation
- Limited to m6A

**Benchmarking Considerations:**
- Feature extraction quality
- SVM model parameters
- Training data match

---

### 9. TandemMod

**Publication:** [Nature Communications 2024](https://www.nature.com/articles/s41467-024-48437-4)

**Input Requirements:**
- Filtered BAM files
- BAM index files
- Transcriptome FASTA
- Configuration YAML

**Output Format:**
```
transcript    position    modification    probability    stoichiometry
```

**Supported Modifications:** Multiple (transfer learning framework)

**Pros:**
- Deep learning with transfer learning
- Multi-modification detection
- Single-sample mode (no control needed)
- Provides stoichiometry estimates

**Cons:**
- Requires GPU for training
- Model selection affects results
- Newer tool with less validation

**Benchmarking Considerations:**
- Model type selection (single/multi)
- Probability threshold tuning
- Transfer learning quality

---

### 10. DirectRM

**Input Requirements:**
- Filtered BAM files
- BAM index files
- Transcriptome FASTA
- Configuration YAML

**Output Format:**
```
transcript    position    modification    probability
```

**Supported Modifications:** Multiple

**Pros:**
- Single-sample detection
- Multiple modification types
- Configurable thresholds

**Cons:**
- Less documented
- Model-dependent
- Requires configuration

**Benchmarking Considerations:**
- min_prob threshold (default 0.7)
- Modification type selection
- Coverage requirements

---

### 11. m6ATM

**Publication:** [Nature Communications 2023](https://www.nature.com/articles/s41467-023-37596-5)

**Input Requirements:**
- Filtered BAM files
- BAM index files
- Transcriptome FASTA
- Configuration YAML

**Output Format:**
```
transcript    position    probability    stoichiometry
```

**Supported Modifications:** m6A

**Pros:**
- Transformer-based architecture
- High performance in benchmarks
- Provides stoichiometry estimates
- Single-sample mode

**Cons:**
- Limited to m6A
- Memory-intensive
- Requires specific threshold tuning

**Benchmarking Considerations:**
- Probability threshold (default 0.6)
- Stoichiometry calculation accuracy
- Coverage depth impact

---

### 12. Rnano

**Input Requirements:**
- Filtered BAM files
- BAM index files
- Transcriptome FASTA
- Configuration YAML

**Output Format:**
```
transcript    position    modification    probability
```

**Supported Modifications:** m6A

**Pros:**
- Pretrained model available
- GPU optional
- Single-sample detection

**Cons:**
- Less documented
- Limited modification types

**Benchmarking Considerations:**
- Model selection
- CPU vs GPU performance
- Threshold tuning

---

### 13. NanoPSU

**GitHub:** [sihaohuanguc/Nanopore_psU](https://github.com/sihaohuanguc/Nanopore_psU)

**Input Requirements:**
- Filtered BAM files
- BAM index files
- Transcriptome FASTA

**Output Format:**
```
transcript    position    psi_probability    coverage
```

**Supported Modifications:** Pseudouridine (Ψ)

**Pros:**
- Specialized for Ψ detection
- No minimum read requirement
- Single-sample mode

**Cons:**
- Limited to pseudouridine
- Newer tool

**Benchmarking Considerations:**
- Ψ-specific signal patterns
- Coverage thresholds
- False positive rates

---

### 14. NanoMUD

**Publication:** [BioRxiv 2024](https://www.biorxiv.org/content/10.1101/2024.05.08.593203)

**Input Requirements:**
- Filtered BAM files
- BAM index files
- Transcriptome FASTA

**Output Format:**
```
transcript    position    modification    probability    coverage
```

**Supported Modifications:** Ψ, m1Ψ (configurable)

**Pros:**
- Multiple Ψ variants
- Configurable modification list
- Single-sample mode

**Cons:**
- Limited to Ψ-related modifications
- Newer tool

**Benchmarking Considerations:**
- Modification type parameter
- Probability thresholds
- Validation against known Ψ sites

---

### 15. Penguin

**Input Requirements:**
- Filtered BAM files
- BAM index files
- Transcriptome FASTA

**Output Format:**
```
transcript    position    psi_probability    score
```

**Supported Modifications:** Pseudouridine (Ψ)

**Pros:**
- Specialized Ψ detection
- Single-sample mode
- Fast execution

**Cons:**
- Limited to pseudouridine
- Less documentation

**Benchmarking Considerations:**
- Comparison with NanoPSU/NanoMUD
- Sensitivity vs specificity trade-offs
- Coverage requirements

---

### 16. psipore

**Input Requirements:**
- Filtered BAM files (Native and Control)
- BAM index files
- Transcriptome FASTA
- Regions BED file

**Output Format:**
```
transcript    position    psi_score    p_value    coverage
```

**Supported Modifications:** Pseudouridine (Ψ)

**Pros:**
- Comparative approach for Ψ
- Region-specific analysis
- Statistical framework

**Cons:**
- Requires control sample
- Region preparation needed
- Limited to Ψ

**Benchmarking Considerations:**
- Comparative vs single-sample Ψ tools
- Region overlap impact
- Statistical significance thresholds

---

## Summary Comparison Table

| Tool | Mode | Input | Signal Source | m6A | Ψ | Other | Memory |
|------|------|-------|---------------|-----|---|-------|--------|
| xPore | Comparative | eventalign | f5c | ✓ | - | - | High |
| Nanocompore | Comparative | eventalign | f5c | ✓ | - | - | High |
| Baleen | Comparative | eventalign | f5c | ✓ | - | - | Very High |
| pyBaleen | Comparative | BAM+BLOW5 | raw | ✓ | - | - | Very High |
| DiffErr | Comparative | BAM | errors | ✓ | - | - | Medium |
| DRUMMER | Comparative | BAM | errors | ✓ | - | ✓ | Low |
| ELIGOS2 | Comparative | BAM | errors | ✓ | ✓ | ✓ | Medium |
| EpiNano | Comparative | BAM (CSV) | errors | ✓ | - | - | Medium |
| psipore | Comparative | BAM | errors | - | ✓ | - | Medium |
| TandemMod | Single-sample | BAM | DL model | ✓ | ✓ | ✓ | High |
| DirectRM | Single-sample | BAM | ML model | ✓ | - | ✓ | High |
| m6ATM | Single-sample | BAM | Transformer | ✓ | - | - | High |
| Rnano | Single-sample | BAM | pretrained | ✓ | - | - | High |
| NanoPSU | Single-sample | BAM | ML model | - | ✓ | - | Medium |
| NanoMUD | Single-sample | BAM | ML model | - | ✓ | ✓ | Medium |
| Penguin | Single-sample | BAM | ML model | - | ✓ | - | Medium |

---

## Benchmarking Recommendations

### Key Metrics to Evaluate

1. **Sensitivity (Recall):** True positive rate for known modification sites
2. **Specificity:** True negative rate for unmodified sites
3. **Precision:** Positive predictive value
4. **F1 Score:** Harmonic mean of precision and recall
5. **MCC (Matthews Correlation Coefficient):** Balanced measure for imbalanced data
6. **AUC-ROC:** Area under ROC curve
7. **Stoichiometry Accuracy:** For tools providing modification frequency estimates

### Benchmarking Dataset Requirements

1. **Synthetic RNA controls:** IVT transcripts with known modifications
2. **Cell line standards:** HEK293T, Hela with validated modification sites
3. **Negative controls:** Unmodified IVT RNA
4. **Coverage variation:** Test at different sequencing depths

### Critical Parameters to Test

| Tool Category | Key Parameters |
|---------------|----------------|
| Signal-based | min_depth, signal normalization, kmer context |
| Error-based | error_rate_threshold, coverage_min |
| ML-based | probability_threshold, model_selection |

### Common Pitfalls

1. **Batch effects:** Different flow cells, chemistry versions
2. **Coverage bias:** Low coverage regions have higher false positives
3. **Sequence context:** Some tools perform poorly at specific kmers
4. **Control matching:** Comparative tools need well-matched controls
5. **Model applicability:** ML tools trained on specific data may not generalize

---

## References

1. [Systematic comparison of tools for m6A mapping from nanopore DRS](https://www.nature.com/articles/s41467-023-37596-5) - Nature Communications 2023
2. [Comparative evaluation of computational models for RNA modification detection](https://academic.oup.com/bib/article/26/4/bbaf404/8233723) - Briefings in Bioinformatics 2024
3. [Nanopore DRS for RNA modification analysis](https://link.springer.com/article/10.1007/s44307-025-00093-5) - 2025 Review
4. [TandemMod: Transferable deep learning framework](https://www.nature.com/articles/s41467-024-48437-4) - Nature Communications 2024
5. [RNA modifications detection by comparative Nanopore DRS](https://www.nature.com/articles/s41467-021-27393-3) - Nature Communications 2021
6. [RMaP Challenge](https://www.nature.com/articles/s42004-025-01507-0) - Nature Communications Chemistry 2025
