# Singularity-Only Mode - Container Configuration

All `conda:` directives have been removed from the workflow. The workflow now uses singularity containers only.

## Container Configuration (config/config.yaml)

| 容器名称 | Docker 镜像 | 用途 | 状态 |
|---------|------------|------|------|
| `python3` | `btrspg/misc-python3:2026-03-20` | Python 3 脚本、数据处理、可视化 | ✅ 就绪 |
| `r_viz` | `6oclock/r-deseq2-edger-ggplot2:3.6.3` | R 绘图 (ggplot2) | ✅ 就绪 |
| `picard` | `broadinstitute/picard:3.1.0` | 序列字典创建 | ✅ 就绪 |

## How to Run

```bash
# Singularity-only mode (no conda)
snakemake --use-singularity --cores 4 -p
```

## Summary of Changes Made

1. **All rule files**: Removed all `conda:` directives (20 files)
2. **prep_reference.smk**: Replaced picard wrapper with direct shell command using `get_container("picard")`
3. **config/config.yaml**: Consolidated pandas/matplotlib containers into python3 (same image, reduces config complexity)
