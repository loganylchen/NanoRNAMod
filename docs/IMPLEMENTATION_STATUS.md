# Implementation Progress

## Phase 1: Multi-Modification Tools ✓ COMPLETE

### TandemMod
- [x] `workflow/envs/tandemmod.yaml` - Conda environment created
- [x] `workflow/rules/modetect_tandemmod.smk` - Rule file created
- [x] `workflow/scripts/tandemmod_config.py` - Config script created
- [x] `workflow/scripts/tandemmod_postprocess.py` - Postprocess script created
- [x] `workflow/scripts/format.py` - Updated with format_tandemmod()
- [x] `workflow/rules/post_format.smk` - Added post_tandemmod and tandemmod_annotate rules
- [x] `workflow/rules/common.smk` - Updated get_final_output()
- [x] `workflow/Snakefile` - Added include
- [x] `config/config.yaml` - Added params, threads, and tools activation

### DirectRM
- [x] `workflow/envs/directrm.yaml` - Conda environment created
- [x] `workflow/rules/modetect_directrm.smk` - Rule file created
- [x] `workflow/scripts/directrm_config.py` - Config script created
- [x] `workflow/scripts/directrm_postprocess.py` - Postprocess script created
- [x] `workflow/scripts/format.py` - Updated with format_directrm()
- [x] `workflow/rules/post_format.smk` - Added post_directrm and directrm_annotate rules
- [x] `workflow/rules/common.smk` - Updated get_final_output()
- [x] `workflow/Snakefile` - Added include
- [x] `config/config.yaml` - Added params, threads, and tools activation

## Phase 2: Enhanced m6A Detection ✓ COMPLETE

### m6ATM
- [x] `workflow/envs/m6atm.yaml` - Conda environment created
- [x] `workflow/rules/modetect_m6atm.smk` - Rule file created
- [x] `workflow/scripts/m6atm_config.py` - Config script created
- [x] `workflow/scripts/m6atm_postprocess.py` - Postprocess script created
- [x] `workflow/scripts/format.py` - Updated with format_m6atm()
- [x] `workflow/rules/post_format.smk` - Added post_m6atm and m6atm_annotate rules
- [x] `workflow/rules/common.smk` - Updated get_final_output()
- [x] `workflow/Snakefile` - Added include
- [x] `config/config.yaml` - Added params, threads, and tools activation

### RNANO
- [x] `workflow/envs/rnano.yaml` - Conda environment created
- [x] `workflow/rules/modetect_rnano.smk` - Rule file created
- [x] `workflow/scripts/rnano_config.py` - Config script created
- [x] `workflow/scripts/rnano_postprocess.py` - Postprocess script created
- [x] `workflow/scripts/format.py` - Updated with format_rnano()
- [x] `workflow/rules/post_format.smk` - Added post_rnano and rnano_annotate rules
- [x] `workflow/rules/common.smk` - Updated get_final_output()
- [x] `workflow/Snakefile` - Added include
- [x] `config/config.yaml` - Added params, threads, and tools activation

## Summary

### Files Created (18 total)
**Environment Files (4)**:
- workflow/envs/tandemmod.yaml
- workflow/envs/directrm.yaml
- workflow/envs/m6atm.yaml
- workflow/envs/rnano.yaml

**Rule Files (4)**:
- workflow/rules/modetect_tandemmod.smk
- workflow/rules/modetect_directrm.smk
- workflow/rules/modetect_m6atm.smk
- workflow/rules/modetect_rnano.smk

**Helper Scripts (8)**:
- workflow/scripts/tandemmod_config.py
- workflow/scripts/tandemmod_postprocess.py
- workflow/scripts/directrm_config.py
- workflow/scripts/directrm_postprocess.py
- workflow/scripts/m6atm_config.py
- workflow/scripts/m6atm_postprocess.py
- workflow/scripts/rnano_config.py
- workflow/scripts/rnano_postprocess.py

**Updated Files (4)**:
- workflow/scripts/format.py (added 2 new format functions)
- workflow/rules/post_format.smk (added 8 new rules)
- workflow/rules/common.smk (added to get_final_output())
- workflow/Snakefile (added 2 includes)
- config/config.yaml (added params, threads, and tools)

## Phase 3: Pseudouridine Detection ⚠️ PENDING

**Tools to Implement (4)**:
1. PsiNanopore - Comparative, highest Ψ accuracy
2. NanoPsu - Standard Ψ tool with stoichiometry
3. NanoMUD - Unique Ψ + N1-methylpseudouridine (m1Ψ) detection
4. Penguin - Ensemble ML models for Ψ

**Estimated Work**:
- 4 conda environment files
- 4 rule files
- 8 helper scripts
- Updates to format.py (4 functions)
- Updates to post_format.smk (8 rules)
- Updates to common.smk (get_final_output)
- Updates to Snakefile (4 includes)
- Updates to config.yaml

## Phase 4: ONT Integration & Additional Tools ⚠️ PENDING

**Tools to Implement (2)**:
1. Dorado - ONT native basecaller with RNA004 support
2. EpiNano2 - High-accuracy m6A detection

**Estimated Work**:
- 2 conda environment files
- 2-3 rule files (Dorado may need 2 rules)
- 4-6 helper scripts
- Updates to format.py (2-3 functions)
- Updates to post_format.smk (4-6 rules)
- Updates to common.smk (get_final_output)
- Updates to Snakefile (2-3 includes)
- Updates to config.yaml

## Current Tool Count

**Before Implementation**: 7 tools
**After Phase 1-2**: 11 tools (4 new added)
**After Phase 3-4 (complete)**: 17 tools (6 more needed)

## Active Tools in config.yaml

**Enabled by default (7)**:
- xpore
- nanocompore
- baleen
- differr
- drummer
- eligos2
- epinano

**Disabled by default (4)**:
- tandemmod
- directrm
- m6atm
- rnano

**To enable a tool**, set `activate: True` in config.yaml:
```yaml
tools:
  tandemmod:
    activate: True
```

## Next Steps

### Immediate (Phase 3)
1. Research PsiNanopore GitHub repository
2. Create psipore.yaml environment
3. Create modetect_psipore.smk rule file
4. Create helper scripts (config, predict, postprocess)
5. Repeat for NanoPsu, NanoMUD, Penguin

### Follow-up (Phase 4)
1. Research Dorado installation and RNA004 support
2. Create dorado.yaml environment
3. Create modetect_dorado.smk rule file
4. Optionally create prep_basecall_dorado.smk for RNA004 basecalling
5. Create EpiNano2 implementation
6. Update all configuration files

## Notes

- All new tools follow the existing workflow architecture
- Single-sample tools work on individual BAM files
- Comparative tools require native/control sample pairs
- All outputs are standardized to TSV format
- All tools can be independently enabled/disabled
- Format and annotation use Baleen.py for GTF annotation

## Testing Status

**Not yet tested** - The implemented code structure is complete but requires:
1. Actual tool installations (pip install, GitHub clones)
2. Test data availability
3. Dry-run validation: `snakemake --dry-run`
4. Integration testing with sample data
5. Benchmark resource usage

## Remaining Work

**Phase 3**: 4 tools, ~3 weeks
**Phase 4**: 2 tools, ~2 weeks
**Testing & Documentation**: ~2 weeks

**Total remaining**: ~7 weeks
