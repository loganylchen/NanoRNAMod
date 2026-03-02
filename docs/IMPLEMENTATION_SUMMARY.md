# Implementation Plan Summary

## Quick Reference

This document provides a quick-reference summary of the full implementation plan in `IMPLEMENTATION_PLAN.md`.

---

## Tools to Add (10 Total)

### Phase 1: Multi-Modification Tools (Weeks 1-3)

| Tool | Type | Mods | Priority | Files to Create |
|-------|-------|-------|-----------------|
| **TandemMod** | Single-sample | High | 1 rule, 1 env, 2 scripts |
| **DirectRM** | Single-sample | High | 1 rule, 1 env, 2 scripts |

**Total Phase 1**: 2 rules, 2 envs, 4 scripts, 3 weeks

### Phase 2: Enhanced m6A Detection (Weeks 4-5)

| Tool | Type | Mods | Priority | Files to Create |
|-------|-------|-------|-----------------|
| **m6ATM** | Single-sample | High | 1 rule, 1 env, 1 script |
| **RNANO** | Single-sample | High | 1 rule, 1 env, 2 scripts |

**Total Phase 2**: 2 rules, 2 envs, 3 scripts, 2 weeks

### Phase 3: Pseudouridine Detection (Weeks 6-8)

| Tool | Type | Mods | Priority | Files to Create |
|-------|-------|-------|-----------------|
| **PsiNanopore** | Comparative | High | 1 rule, 1 env, 1 script |
| **NanoPsu** | Single-sample | High | 1 rule, 1 env, 1 script |
| **NanoMUD** | Single-sample | High | 1 rule, 1 env, 1 script |
| **Penguin** | Single-sample | High | 1 rule, 1 env, 1 script |

**Total Phase 3**: 4 rules, 4 envs, 4 scripts, 3 weeks

### Phase 4: ONT Integration (Weeks 9-10)

| Tool | Type | Mods | Priority | Files to Create |
|-------|-------|-------|-----------------|
| **Dorado** | Single-sample (basecalling) | Medium | 2 rules, 1 env, 2 scripts |
| **EpiNano2** | Comparative | Medium | 1 rule, 1 env, 1 script |

**Total Phase 4**: 3 rules, 2 envs, 4 scripts, 2 weeks

---

## Files to Create Summary

### Rule Files (10)
```
workflow/rules/modetect_tandemmod.smk
workflow/rules/modetect_directrm.smk
workflow/rules/modetect_m6atm.smk
workflow/rules/modetect_rnano.smk
workflow/rules/modetect_psipore.smk
workflow/rules/modetect_nanopsu.smk
workflow/rules/modetect_nanomud.smk
workflow/rules/modetect_penguin.smk
workflow/rules/modetect_dorado.smk
workflow/rules/modetect_epinano2.smk
```

**Optional:**
```
workflow/rules/prep_basecall_dorado.smk  # If using Dorado for basecalling
```

### Conda Environment Files (10)
```
workflow/envs/tandemmod.yaml
workflow/envs/directrm.yaml
workflow/envs/m6atm.yaml
workflow/envs/rnano.yaml
workflow/envs/psipore.yaml
workflow/envs/nanopsu.yaml
workflow/envs/nanomud.yaml
workflow/envs/penguin.yaml
workflow/envs/dorado.yaml
workflow/envs/epinano2.yaml
```

### Helper Scripts (15-20)
```
workflow/scripts/tandemmod_config.py
workflow/scripts/tandemmod_predict.py
workflow/scripts/directrm_config.py
workflow/scripts/directrm_predict.py
workflow/scripts/m6atm_predict.py
workflow/scripts/rnano_predict.py
workflow/scripts/psipore_config.py
workflow/scripts/nanopsu_predict.py
workflow/scripts/nanomud_predict.py
workflow/scripts/penguin_predict.py
workflow/scripts/dorado_download_models.py
workflow/scripts/dorado_predict.py
workflow/scripts/epinano2_config.py
workflow/scripts/download_tools.py
```

---

## Configuration Updates

### config.yaml Additions

```yaml
params:
  # New parameters (add to existing)
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
  # New threads (add to existing)
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
  # New tools (add to existing)
  tandemmod:
    activate: False
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
    basecaller: False
  epinano2:
    activate: False
```

### workflow/Snakefile Additions

```python
# Add these includes after existing includes
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

---

## Testing Checklist

### Before Starting Each Phase

- [ ] Read tool documentation thoroughly
- [ ] Locate GitHub repository
- [ ] Identify dependencies
- [ ] Test installation locally
- [ ] Create conda environment file
- [ ] Create helper scripts
- [ ] Create rule file
- [ ] Test rule with dry-run
- [ ] Test rule with small dataset

### During Implementation

- [ ] Follow project coding conventions (AGENTS.md)
- [ ] Use KEEP_OR_NOT for intermediate files
- [ ] Define log files
- [ ] Define benchmark files
- [ ] Set resources (memory, threads)
- [ ] Set priority
- [ ] Use conda environment
- [ ] Redirect stderr to log
- [ ] No inline comments (project convention)

### After Each Tool

- [ ] Rule executes successfully
- [ ] Output file exists
- [ ] Output format is correct
- [ ] Required columns present
- [ ] Integration with post_format works
- [ ] Annotation works (if GTF provided)
- [ ] Resource usage documented
- [ ] Documentation updated

---

## Key Implementation Notes

### Comparative vs Single-Sample

**Comparative Tools** (require native/control samples):
- PsiNanopore (NEW)
- EpiNano2 (NEW)
- xpore (existing)
- baleen (existing)
- nanocompore (existing)
- differr (existing)
- drummer (existing)
- eligos2 (existing)
- epinano (existing)

**Single-Sample Tools** (work on individual samples):
- TandemMod (NEW)
- DirectRM (NEW)
- m6ATM (NEW)
- RNANO (NEW)
- NanoPsu (NEW)
- NanoMUD (NEW)
- Penguin (NEW)

### Input Requirements

**Eventalign-Dependent** (require f5c eventalign):
- baleen (existing)
- nanocompore (existing)

**Raw Signal-Based** (no eventalign needed):
- TandemMod (NEW)
- DirectRM (NEW)
- m6ATM (NEW)
- RNANO (NEW)
- Dorado (NEW - basecalling)
- NanoPsu (NEW)
- NanoMUD (NEW)
- Penguin (NEW)
- PsiNanopore (NEW - uses BAM)
- EpiNano2 (NEW - uses BAM)

### Special Features

**Unique Capabilities**:
- **TandemMod**: Transfer learning, 3 mods
- **DirectRM**: Most comprehensive (6 mods), crosstalk detection
- **RNANO**: NLP approach, transformers
- **m6ATM**: Better than m6Anet, stoichiometry
- **NanoMUD**: Only tool detecting m1Ψ
- **Dorado**: ONT native, RNA004 chemistry support, real-time
- **PsiNanopore**: Highest Ψ accuracy, hypermodifications

---

## Common Patterns

### Rule File Template

```smk
rule {tool}_{action}:
    input:
        # Tool-specific inputs
    output:
        # Standardized output
    params:
        # Configuration parameters
    threads: config["threads"]["{tool}"]
    resources:
        mem_mb = <memory_requirement>
    priority: 10
    log:
        "logs/{project}/{tool}/{wildcards}.log"
    conda:
        "../envs/{tool}.yaml"
    shell/script:
        # Tool execution
```

### Helper Script Template

```python
import sys
import os
import pandas as pd

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = open(snakemake.log[1], "w")

# Snakemake parameters
tool_param = snakemake.params.tool_param
input_file = snakemake.input.input_file
output_dir = snakemake.output.output_dir

# Tool-specific logic
df = pd.read_csv(input_file, sep='\t')
# Process data...
df.to_csv(output_file, sep='\t', index=False)
```

### Format Function Template

```python
def format_{tool}(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    df = df.sort_values(['transcript', 'position'])
    # Tool-specific formatting
    df.to_csv(output_file, sep='\t', index=False)
```

---

## Quick Commands

### Testing Single Tool

```bash
# Enable tool in config.yaml
# tools:
#   tandemmod:
#     activate: True

# Run specific tool
cd .test
snakemake --use-conda --cores 2 -R tandemmod_predict
snakemake --use-conda --cores 2 -R post_tandemmod
snakemake --use-conda --cores 2 -R tandemmod_annotate
```

### Testing Multiple Tools

```bash
# Enable multiple tools in config.yaml

# Run full workflow
snakemake --use-conda --cores 4 -k
```

### Dry Run

```bash
# Check dependencies without executing
snakemake --dry-run --cores 1
```

---

## Progress Tracking

### Phase 1 (Weeks 1-3)
- [ ] TandemMod implemented
- [ ] DirectRM implemented
- [ ] Both tested
- [ ] Both documented

### Phase 2 (Weeks 4-5)
- [ ] m6ATM implemented
- [ ] RNANO implemented
- [ ] Both tested
- [ ] Both documented

### Phase 3 (Weeks 6-8)
- [ ] PsiNanopore implemented
- [ ] NanoPsu implemented
- [ ] NanoMUD implemented
- [ ] Penguin implemented
- [ ] All tested
- [ ] All documented

### Phase 4 (Weeks 9-10)
- [ ] Dorado implemented
- [ ] EpiNano2 implemented
- [ ] Both tested
- [ ] Both documented

### Final Review (Weeks 11-12)
- [ ] All 10 tools operational
- [ ] Full workflow tested
- [ ] Documentation complete
- [ ] Release ready

---

## Resources and Estimates

### Development Effort

| Phase | Estimated Hours | Complexity |
|-------|-----------------|------------|
| Phase 1 | 80-120 | High |
| Phase 2 | 60-80 | Medium |
| Phase 3 | 80-100 | High |
| Phase 4 | 60-80 | Medium |
| Testing & Documentation | 40-60 | Low |
| **Total** | **320-440 hours** | **8-11 weeks** |

### Code Volume

| Type | Count | Estimated Lines |
|-------|--------|----------------|
| Rule files | 10 | 2,000-2,500 |
| Environment files | 10 | 100-150 |
| Helper scripts | 15-20 | 800-1,200 |
| Updates to existing files | 4 | 200-300 |
| **Total** | **39-44 files** | **3,100-4,150 lines** |

### Testing Effort

| Test Type | Count | Estimated Hours |
|------------|--------|-----------------|
| Unit tests | 10 tools | 20-30 |
| Integration tests | 4 phases | 20-30 |
| Benchmark tests | 10 tools | 15-20 |
| **Total** | - | **55-80 hours** |

---

## Success Metrics

### Quantitative Targets

- [ ] 10/10 tools successfully implemented
- [ ] All tools pass unit tests
- [ ] All tools pass integration tests
- [ ] All tools produce standard output format
- [ ] Documentation is complete
- [ ] Timeline met (10-12 weeks)

### Qualitative Targets

- [ ] Code follows project conventions
- [ ] No dependency conflicts
- [ ] Resource usage is optimized
- [ ] User experience is smooth
- [ ] Tool comparison is easy

---

## Contact and Support

### Tool Authors

Reach out to tool authors if needed:
- Check GitHub Issues
- Read documentation thoroughly
- Check tool-specific publications

### Community Support

- Join tool-specific mailing lists/forums
- Monitor Snakemake community
- Stay updated with tool releases

---

## Next Steps

1. **Review this summary** with team
2. **Approve implementation plan**
3. **Start Phase 1 implementation**
4. **Track progress using checklists**
5. **Adjust timeline as needed**
6. **Document lessons learned**
7. **Prepare for release**

---

**For detailed information, see**: `IMPLEMENTATION_PLAN.md`
**For literature review, see**: `LITERATURE_REVIEW.md`
**For tool comparison, see**: `TOOLS_SUMMARY.md`
