# Nanopore RNA Modification Tools: Summary Table

## Current Tools in NanoRNAMod

| Tool | Year | Type | Modifications | Citations | Status |
|-------|-------|-------|--------------|--------|
| xpore | 2021 | Comparative, multiple | N/A | ✓ Included |
| baleen | 2021 | UMAP+GMM, multiple | N/A | ✓ Included |
| nanocompore | 2021 | Model-free comparative | N/A | ✓ Included |
| differr | 2021 | Statistical differential | N/A | ✓ Included |
| drummer | 2021 | Comparative | N/A | ✓ Included |
| eligos2 | 2021 | Machine learning | N/A | ✓ Included |
| epinano | 2021 | Neural network | N/A | ✓ Included |

## Recommended Additions by Priority

### HIGH PRIORITY (Immediate Implementation)

| Tool | Year | Type | Modifications | Citations | Chemistry | Notes |
|-------|-------|-------|--------------|------------|--------|
| **TandemMod** | 2024 | Transfer learning, multi-mod (m6A, Ψ, m5C) | 67 | RNA002, RNA004 | State-of-the-art, Nature paper |
| **Dorado** | 2024 | Built-in basecaller (m6A, m5C, Ψ) | N/A | RNA004 (native) | Official ONT, real-time |
| **DirectRM** | 2025 | Multi-mod (m6A, Ψ, m5C, A-to-I, m7G, m1A) | 3 | RNA002 | Most comprehensive (6 types) |
| **m6ATM** | 2024 | Deep learning, m6A | 11 | RNA002 | Better than m6Anet |
| **RNANO** | 2025 | NLP, multi-mod | 1 | RNA002 | Novel approach, recent |

### MEDIUM PRIORITY (Short-term Implementation)

| Tool | Year | Type | Modifications | Citations | Chemistry | Notes |
|-------|-------|-------|--------------|------------|--------|
| **PsiNanopore** | 2023 | Comparative, Ψ, hypermodifications | 120 | RNA002 | Highest Ψ accuracy, most cited |
| **NanoPsu** | 2021 | Quantitative, Ψ | 118 | RNA002 | Standard Ψ tool, highly cited |
| **NanoMUD** | 2024 | Regression, Ψ, N1-methylpseudouridine (m1Ψ) | 8 | RNA002 | Only tool for m1Ψ, AUC 0.998 |
| **Penguin** | 2022 | ML, Ψ | 57 | RNA002 | Popular Ψ tool |
| **EpiNano** | 2021 | m6A | 87 | RNA002 | High accuracy (~90%) |

### LOW PRIORITY (Long-term Consideration)

| Tool | Year | Type | Modifications | Citations | Chemistry | Notes |
|-------|-------|-------|--------------|------------|--------|
| **NanoListener** | 2026 | Transformer, multi-mod | N/A | N/A | Emerging tech, audio-based |
| **NanoSpeech** | 2026 | Audio-to-text, multi-mod | N/A | N/A | Emerging tech, early stage |
| **MoDorado** | 2025 | Enhanced tRNA modifications | 9 | RNA004 | tRNA-specific |
| **DeepRM** | 2025 | Deep learning, multi-mod | 2 | RNA002 | Powerful but less established |

## Modification Type Coverage

| Modification | Current Tools | After Adding High Priority |
|-------------|---------------|-------------------------|
| m6A | xpore, baleen, nanocompore, differr, drummer, eligos2, epinano | +TandemMod, Dorado, DirectRM, m6ATM, RNANO, EpiNano |
| Ψ | epinano | +TandemMod, Dorado, DirectRM, RNANO, PsiNanopore, NanoPsu, NanoMUD, Penguin |
| m5C | xpore, baleen, nanocompore | +TandemMod, Dorado, DirectRM, RNANO |
| A-to-I | xpore | +DirectRM, RNANO |
| m7G | xpore | +DirectRM, RNANO |
| m1A | xpore | +DirectRM, RNANO |
| m1Ψ | None | +NanoMUD |

## Tool Characteristics Comparison

| Tool | Comparative? | Single-sample? | Eventalign Required? | Raw Signal? | Stoichiometry? |
|-------|---------------|-----------------|---------------------|-------------|----------------|
| **Current Tools** |
| xpore | ✓ | ✗ | ✗ | ✓ | ✗ |
| baleen | ✓ | ✗ | ✓ | ✗ | ✓ |
| nanocompore | ✓ | ✗ | ✓ | ✗ | ✗ |
| differr | ✓ | ✗ | ✗ | ✗ | ✗ |
| drummer | ✓ | ✗ | ✗ | ✗ | ✗ |
| eligos2 | ✓ | ✗ | ✗ | ✗ | ✗ |
| epinano | ✓ | ✗ | ✗ | ✗ | ✗ |
| **Recommended High Priority** |
| TandemMod | ✗ | ✓ | ✗ | ✓ | ✗ |
| Dorado | ✗ | ✓ | ✗ | ✓ | ✓ |
| DirectRM | ✗ | ✓ | ✗ | ✓ | ✓ |
| m6ATM | ✗ | ✓ | ✗ | ✓ | ✓ |
| RNANO | ✗ | ✓ | ✗ | ✓ | ✗ |
| **Recommended Medium Priority** |
| PsiNanopore | ✓ | ✗ | ✗ | ✗ | ✓ |
| NanoPsu | ✗ | ✓ | ✗ | ✓ | ✓ |
| NanoMUD | ✗ | ✓ | ✗ | ✓ | ✓ |
| Penguin | ✗ | ✓ | ✗ | ✓ | ✗ |
| EpiNano | ✗ | ✓ | ✗ | ✓ | ✗ |

## Key Insights

### Chemistry Support
- **RNA002**: Most current tools optimized for this chemistry
- **RNA004**: New chemistry (2024) requires RNA-specific flowcell; Dorado native support
- **Future**: RNA004 will become standard; RNA002 tools may need retraining

### Comparative vs Single-Sample
- **Comparative (7 current + 1 recommended)**: Require control samples, generally higher specificity
- **Single-sample (5 high priority)**: No control needed, faster, more flexible

### Novel Approaches
- **NLP-based**: RNANO treats signals as natural language
- **Audio-based**: NanoListener/NanoSpeech treat signals as audio
- **Transfer learning**: TandemMod leverages knowledge from other domains

### Coverage Gaps
- **m1Ψ**: Only NanoMUD detects this modification
- **tRNA modifications**: Only MoDorado specializes in tRNAs
- **Hypermodifications**: PsiNanopore detects type I/II hypermodifications

## Implementation Recommendations

1. **Add all 5 high-priority tools** (TandemMod, Dorado, DirectRM, m6ATM, RNANO)
   - Immediate impact: multi-modification support, RNA004 compatibility
   - Total additions: 5 tools

2. **Add all 5 medium-priority tools** (PsiNanopore, NanoPsu, NanoMUD, Penguin, EpiNano)
   - Enhanced Ψ detection capabilities
   - Coverage of m1Ψ modification
   - Total additions: 5 tools

3. **Evaluate low-priority tools** based on community adoption and validation
   - Monitor citation growth
   - Assess tool stability and maintenance
   - Consider specialized use cases (tRNA, audio-based approaches)

**Total recommended additions: 10 tools**
**Final tool count: 17 tools**
**Modification types covered: 8+**
