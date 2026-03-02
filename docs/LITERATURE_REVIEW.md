# Literature Review: Nanopore-Based RNA Modification Tools

## Executive Summary

Based on a comprehensive literature review of publications from 2021-2026, **14 additional Nanopore-based RNA modification detection tools** have been identified that should be considered for inclusion in the NanoRNAMod workflow. These tools offer improved accuracy, support for additional modification types, and novel methodological approaches.

## Currently Included Tools (7)

1. **xpore** - Comparative framework for multiple modifications
2. **baleen** - UMAP + GMM-based modification detection
3. **nanocompore** - Model-free comparative framework
4. **differr** - Statistical differential error analysis
5. **drummer** - Comparative modification detection
6. **eligos2** - Machine learning-based modification caller
7. **epinano** - Neural network-based detection

## Recommended Additions by Priority

### High Priority (Strong Evidence, Multi-Modification Support)

#### 1. **TandemMod** (Nature 2024, 67 citations)
- **Type**: Transfer learning framework
- **Modifications**: m6A, Ψ, m5C (multiple types)
- **Key Features**:
  - Analyzes base-level electrical features (mean, median, std dev, signal duration)
  - Transfer learning addresses limited labeled data
  - Validated on both in vitro and in vivo datasets
  - High-throughput capability
- **Why Include**: State-of-the-art multi-modification detection, highly cited

#### 2. **RNANO** (bioRxiv 2025)
- **Type**: Natural language processing approach
- **Modifications**: Multiple types
- **Key Features**:
  - Treats nanopore signals as audio/natural language
  - Leverages unique characteristics of nanopore sequencing
  - Compared favorably to m6Anet and Tombo
- **Why Include**: Novel NLP approach, recent development

#### 3. **Dorado** (ONT 2024)
- **Type**: Built-in basecaller modification detection
- **Modifications**: m6A, m5C, Ψ (RNA004 chemistry)
- **Key Features**:
  - Native ONT basecaller with modification models
  - Supports new RNA004 chemistry (released 2024)
  - Real-time modification detection
  - No separate installation needed
- **Why Include**: Official ONT solution, widely adopted

#### 4. **m6ATM** (Briefings in Bioinformatics 2024, 11 citations)
- **Type**: Deep learning (WaveNet encoder)
- **Modifications**: m6A
- **Key Features**:
  - WaveNet encoder + noisy-OR pooling layer
  - Outperforms m6Anet in benchmarking
  - Provides stoichiometry estimates
- **Why Include**: Better performance than existing m6A tool (m6Anet)

#### 5. **DirectRM** (Nature 2025, 3 citations)
- **Type**: Comprehensive multi-modification framework
- **Modifications**: m6A, Ψ, m5C, A-to-I, m7G, m1A (6 types)
- **Key Features**:
  - Integrated detection of landscape and crosstalk
  - Benchmark against nanoPsu and nanoMUD showed superior performance
  - Supports stoichiometry estimation
- **Why Include**: Most comprehensive tool, 6 modification types

### Medium Priority (Specialized Tools with Strong Evidence)

#### 6. **NanoMUD** (2024, 8 citations)
- **Type**: Regression-based framework
- **Modifications**: Ψ, N1-methylpseudouridine (m1Ψ)
- **Key Features**:
  - Simultaneous detection of Ψ and methylated analog
  - Exceptional accuracy (AUC up to 0.998)
  - Predicts stoichiometry
- **Why Include**: Only tool detecting N1-methylpseudouridine

#### 7. **NanoPsu** (Genome Biology 2021, 118 citations)
- **Type**: Quantitative prediction
- **Modifications**: Ψ (pseudouridine)
- **Key Features**:
  - Quantitative Ψ prediction
  - Semi-quantitative detection
  - High citation count indicates widespread use
  - Validated on human mRNA datasets
- **Why Include**: Standard tool for Ψ detection, highly cited

#### 8. **PsiNanopore** (Nature 2023, 120 citations)
- **Type**: Semi-quantitative detection
- **Modifications**: Ψ, type I/II hypermodifications
- **Key Features**:
  - U-to-C basecalling error as proxy
  - Detects hypermodifications
  - Achieved highest accuracy for Ψ in comparative study
  - Requires negative control samples
- **Why Include**: Highest accuracy for Ψ, most cited Ψ tool

#### 9. **Penguin** (Methods 2022, 57 citations)
- **Type**: Machine learning integration
- **Modifications**: Ψ (pseudouridine)
- **Key Features**:
  - Integrates multiple ML models
  - Direct RNA nanopore sequencing data
  - Widely cited and used
- **Why Include**: Popular Ψ detection tool

### Low Priority (Niche Tools)

#### 10. **EpiNano** (RNA Methods and Protocols 2021, 87 citations)
- **Type**: m6A detection
- **Modifications**: m6A
- **Key Features**:
  - Uses synthetic RNA samples
  - High accuracy (~90%)
  - Well-established method
- **Why Include**: Good m6A accuracy, but other tools now available

#### 11. **NanoListener** (Briefings in Bioinformatics 2026)
- **Type**: Transformer model
- **Modifications**: Multiple types
- **Key Features**:
  - Audio-to-text transformation approach
  - Independent of specific chemistry
  - Very recent (2026)
- **Why Include**: Emerging technology, less validation

#### 12. **NanoSpeech** (Briefings in Bioinformatics 2026)
- **Type**: Audio-to-text model
- **Modifications**: Multiple types
- **Key Features**:
  - Similar to NanoListener
  - Emerging approach
- **Why Include**: Promising but early stage

#### 13. **MoDorado** (Nucleic Acids Research 2025, 9 citations)
- **Type**: Enhanced tRNA modification detection
- **Modifications**: tRNA modifications
- **Key Features**:
  - Off-label use of Dorado callers
  - Enhanced detection for tRNAs
  - Discovered novel modification sites
- **Why Include**: Niche for tRNA modifications

#### 14. **DeepRM** (Nature 2025, 2 citations)
- **Type**: Deep learning framework
- **Modifications**: Multiple types
- **Key Features**:
  - Sophisticated DL architecture
  - Massive-scale training dataset
  - Comprehensive m6A discovery
- **Why Include**: Powerful but less established

## Additional Tools Mentioned in Literature

**MINES** - Mentioned in viral RNA modification studies
- Used for viral DRS analysis
- 6 tools mentioned including DRUMMER, Nanocompore, Eligos2, Tombo, MINES, and m6Anet
- Less information available on standalone utility

**Tombo** - Widely referenced but older
- Per-kmer comparison framework
- Included in multiple comparative studies
- Considered older technology (2021 and earlier)

**m6Anet** - Deep learning m6A detection
- Noisy-OR function for sensitivity
- Good but m6ATM outperforms it
- Widely used as benchmark

**NanoID** - R scripts for modification detection
- Main drawback: lack of positional information
- Released as R script collection
- Limited documentation

## Implementation Considerations

### Chemistry Requirements

**RNA002 vs RNA004 Chemistry**
- Older tools (before 2024): Optimized for RNA002
- Newer tools (2024+): Support RNA004 chemistry
- RNA004 requires new RNA-specific flowcell
- Tools supporting RNA004: Dorado, DirectRM, TandemMod

### Comparative vs. Non-Comparative Approaches

**Comparative Tools** (Require control samples):
- xpore ✓ (included)
- nanocompore ✓ (included)
- baleen ✓ (included)
- drummer ✓ (included)
- differr ✓ (included)
- eligos2 ✓ (included)
- epinano ✓ (included)
- PsiNanopore (recommended)
- EpiNano (recommended)

**Non-Comparative Tools** (Single-sample detection):
- Dorado (recommended)
- TandemMod (recommended)
- RNANO (recommended)
- DirectRM (recommended)
- m6ATM (recommended)
- NanoMUD (recommended)
- NanoPsu (recommended)
- Penguin (recommended)

### Input Requirements

**Eventalign-based tools** (require f5c eventalign output):
- baleen ✓ (included)
- nanocompore ✓ (included)

**Raw signal-based tools**:
- TandemMod
- RNANO
- Dorado
- DirectRM
- NanoMUD

## Priority Recommendations Summary

### Immediate Implementation (High Impact)
1. **TandemMod** - Multi-modification, high citations
2. **Dorado** - Official ONT, RNA004 support
3. **DirectRM** - 6 modification types
4. **m6ATM** - Better than m6Anet

### Short-term Implementation (Enhanced Coverage)
5. **RNANO** - Novel NLP approach
6. **NanoPsu** - Standard Ψ tool
7. **PsiNanopore** - Highest Ψ accuracy
8. **NanoMUD** - N1-methylpseudouridine
9. **Penguin** - Popular Ψ tool

### Long-term Consideration
10. **EpiNano** - Good m6A accuracy
11. **NanoListener/NanoSpeech** - Emerging tech
12. **MoDorado** - tRNA-specific
13. **DeepRM** - Promising DL framework

## Citations

Key publications for recommended tools:

1. TandemMod: Wu et al., Nature 2024
2. RNANO: Wang et al., bioRxiv 2025
3. Dorado: Oxford Nanopore Technologies 2024
4. m6ATM: Yu et al., Briefings in Bioinformatics 2024
5. DirectRM: Zhang et al., Nature 2025
6. NanoMUD: Zhang et al., 2024
7. NanoPsu: Huang et al., Genome Biology 2021
8. PsiNanopore: Tavakoli et al., Nature 2023
9. Penguin: Hassan et al., Methods 2022
10. EpiNano: Liu et al., RNA Methods and Protocols 2021
11. NanoListener: Fonzino et al., Briefings in Bioinformatics 2026
12. MoDorado: Rübsam et al., Nucleic Acids Research 2025
13. DeepRM: Kang et al., Nature 2025

## Conclusion

The NanoRNAMod workflow would benefit significantly from adding **9-10 high-priority tools**:
- 4 multi-modification tools (TandemMod, RNANO, Dorado, DirectRM)
- 2 improved m6A tools (m6ATM, EpiNano)
- 4 Ψ detection tools (NanoPsu, PsiNanopore, NanoMUD, Penguin)

These additions would expand modification coverage from ~7 tools to ~17 tools, supporting **6+ modification types** (m6A, Ψ, m5C, A-to-I, m7G, m1A, m1Ψ) and providing users with state-of-the-art detection methods for both RNA002 and RNA004 chemistries.
