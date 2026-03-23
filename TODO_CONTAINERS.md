# Singularity-Only Mode - Container Configuration

All `conda:` directives have been removed from the workflow. The workflow now uses singularity containers only.

## Container Configuration (config/config.yaml)

| 容器名称 | Docker 镜像 | 用途 | 状态 |
|---------|------------|------|------|
| `python3` | `btrspg/misc-python3:2026-03-20` | Python 3 脚本 | ✅ 就绪 |
| `pandas` | `btrspg/misc-python3:2026-03-20` | 数据处理 (pandas, numpy) | ✅ 就绪 |
| `matplotlib` | `btrspg/misc-python3:2026-03-20` | 可视化 (matplotlib, seaborn) | ✅ 就绪 |
| `r_viz` | `6oclock/r-deseq2-edger-ggplot2:3.6.3` | R 绘图 (ggplot2) | ✅ 就绪 |
| `picard` | `broadinstitute/picard:3.1.0` | 序列字典创建 | ✅ 就绪 |

## How to Run

```bash
# Singularity-only mode (no conda)
snakemake --use-singularity --cores 4 -p
```

## Summary of Changes Made

1. **All rule files**: Removed all `conda:` directives (20 files)
2. **benchmark_viz.smk**: Changed `container: None` to `container: get_container("tool_name")`
3. **benchmark_accuracy.smk**: Added `container: get_container("pandas")`
4. **benchmark_report.smk**: Added `container: get_container("pandas")`
5. **prep_reference.smk**: Replaced picard wrapper with direct shell command using `get_container("picard")`
6. **config/config.yaml**: Updated container definitions with user-specified images
