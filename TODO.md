# TODO — NanoRNAMod Improvements

## Critical

### 1. samtools 容器指向了 minimap2 镜像
- **File**: `config/config.yaml` line 103
- **Problem**: `samtools: "docker://btrspg/minimap2:2.28"` — samtools 规则实际运行的是 minimap2 容器。虽然该容器恰好包含 samtools，但这容易造成混淆且不可靠。
- **Fix**: 替换为专用 samtools 容器（如 `docker://biocontainers/samtools:1.20` 或构建专用镜像）

### 2. 死代码规则：`*_epi` 规则无下游消费者
- **Files**: `prep_align_minimap2.smk` lines 55-78, `prep_readsfiltering.smk` lines 55-71
- **Problem**: `minimap2_transcriptome_align_epi` 产生 `{sample}_3.2.4.bam`，`samtools_filter_mapped_epi` 产生 `{sample}_of.bam`。这两个文件不被任何下游规则使用，是旧版 EpiNano 工作流的残留。
- **Fix**: 删除这两个规则，同时删除 `prep_data_prep.smk` 中注释掉的 `link_fastq_old` 规则（lines 20-31）

### 3. pybaleen 使用未过滤的 BAM
- **File**: `modetect_pybaleen.smk` lines 6-7
- **Problem**: 输入是 `{native}.bam` 和 `{control}.bam`（未过滤），其他所有 modetect 工具都使用 `{sample}_filtered.bam`。
- **Fix**: 改为 `{native}_filtered.bam` 和 `{control}_filtered.bam`

---

## High Priority

### 4. `format.py` 缺少错误处理
- **File**: `workflow/scripts/format.py`
- **Problem**: 作为所有工具输出格式化的核心脚本，没有任何 try/except。如果工具输出文件格式异常（空文件、缺少列、非数值数据），整个 pipeline 会崩溃。且没有 log 重定向。
- **Fix**:
  - 添加 `sys.stderr = open(snakemake.log[0], "w")` 或等效 log 重定向
  - 为每个 `format_*()` 函数添加 try/except
  - 添加空文件/空 DataFrame 检查
  - 对 `df['pos'].astype(int)-1` 做数值类型验证

### 5. 多个后处理脚本缺少 log 重定向
- **Files**: `nanomud_postprocess.py`, `nanopsu_postprocess.py`, `penguin_postprocess.py`, `psipore_postprocess.py`
- **Problem**: 这些脚本完全没有 `snakemake.log` 处理，stdout/stderr 直接丢失。出错时无法排查。
- **Fix**: 在脚本开头添加 `sys.stdout`/`sys.stderr` 重定向到 snakemake log 文件

### 6. 脚本 params 定义但未传递
- **Problem**: 以下规则定义了 `params.extra` 但脚本中未使用：
  - `epinano_prep` → 已修复（`epinano_prep.sh` 现在使用 `${extra}`）
  - `epinano` rule 的 `params.prefix` → `epinano_differr.sh` 已使用
  - 其他规则需逐一排查
- **Fix**: 审查所有规则的 params 是否被对应脚本实际消费，移除未使用的 params 或在脚本中使用

### 7. Benchmark 脚本大量代码重复
- **Files**: `benchmark_*.py`（15+ 文件）
- **Problem**: `tool_from_path()`, `normalize_columns()`, `detect_score_column()`, `TOOL_SCORE_COLUMNS` 等函数/常量在 10+ 文件中重复定义。任何一处修改都需要同步所有文件。
- **Fix**: 将公共函数和常量集中到 `benchmark_utils.py`，其他脚本 `from benchmark_utils import ...`

---

## Medium Priority

### 8. QC 规则使用未过滤的 BAM
- **Files**:
  - `qc_quantification_nanocount.smk` line 3: `{sample}.bam`
  - `qc_variants_bcftools.smk` line 3: `{sample}.bam`
  - `polya_estimate_nanopolish.smk` lines 3, 36: `{sample}.bam` 和 `{sample}.splice.bam`
- **Problem**: 修改检测工具都用 `_filtered.bam`，但 QC 规则用未过滤的。这可能是刻意的（NanoCount/BCFtools 有自己的过滤逻辑），但不一致。
- **Fix**: 确认每个工具是否需要预过滤，如果需要则改为 `_filtered.bam`，并添加注释说明原因

### 9. Config 中混合中英文注释
- **File**: `config/config.yaml` lines 37-39, 66-80
- **Problem**: threads 和 containers 部分有中文注释，其余全英文。
- **Fix**: 统一为英文

### 10. 测试配置不完整
- **File**: `.test/config/config.yaml`
- **Problems**:
  - 缺少 `threads.default` 配置 — 如果工具回退到 default 线程可能报错
  - 缺少 `containers` 部分
  - 9 个工具（tandemmod, directrm, m6atm, rnano, psipore, nanopsu, nanomud, penguin, pybaleen）在主 config 激活但未测试
- **Fix**: 补全 .test config，至少添加所有激活工具的基本配置

### 11. 测试数据缺少 BLOW5 文件
- **Dir**: `.test/data/`
- **Problem**: 只有 fast5，workflow 期望 `blow5/nanopore.blow5`。需要 f5c eventalign 的工具（xpore, nanocompore, baleen）在测试中可能无法运行完整流程。
- **Fix**: 添加测试用 BLOW5 文件，或在测试 config 中禁用需要 BLOW5 的工具

### 12. `post_drummer` 输入使用 glob 而非具体文件
- **File**: `post_format.smk` lines 503-504
- **Problem**: 输入是整个目录，`format_drummer()` 用 `glob.glob()` 找 `summary.txt`。如果目录中有多个文件或文件名变化，会出问题。
- **Fix**: 将 drummer 的实际输出文件名作为具体输入

### 13. 后处理规则的冗余目录输入
- **File**: `post_format.smk`
- **Problem**: `post_epinano`, `post_eligos2`, `post_nanocompore` 等规则同时声明了具体文件和其父目录作为输入。目录输入是多余的。
- **Fix**: 删除冗余的目录输入

### 14. Schema 定义了 config 中不存在的参数
- **File**: `workflow/schemas/config.schema.yaml`
- **Problem**: Schema 定义了 `n_bootstrap`, `alpha`, `fdr_method`, `coverage_bins` 等 benchmark 参数，但实际 config.yaml 中没有。benchmark 脚本使用 `config.get()` 带默认值处理，所以不会崩溃，但 schema 与实际不一致。
- **Fix**: 在 config.yaml 中添加这些参数（带默认值），或从 schema 中移除

---

## Low Priority

### 15. 日志命名模式不统一
- **Problem**:
  - 有的用 `log: stdout=..., stderr=...`（xpore_run, nanocompore）
  - 有的用 `log: out=..., err=...`（baleen_*）
  - 有的用单文件 `log: "path.log"`（大部分规则）
- **Fix**: 统一为一种模式（建议：单文件 log，脚本内部自行分离 stdout/stderr）

### 16. `epinano_slide_variants.py` 是空文件
- **File**: `workflow/scripts/epinano_slide_variants.py`
- **Problem**: 只有 1 行，无任何功能。
- **Fix**: 删除，或在需要时实现滑动窗口功能

### 17. Legacy conda 环境文件
- **Files**: `workflow/envs/tombo.yaml`, `workflow/envs/m6anet.yaml`, `workflow/envs/ont-fast5-api.yaml`
- **Problem**: 存在但无对应规则，也无对应工具配置。
- **Fix**: 删除或在 config 中添加对应工具

### 18. README 引用不存在的文件
- **File**: `README.md`
- **Problem**: 安装说明引用 `workflow/envs/main.yaml` 但该文件不存在。关于 "uncomment tools in Snakefile" 的说明已过时（所有工具已默认 include）。
- **Fix**: 更新安装说明，移除过时内容

### 19. `samples.tsv` 的 Directory 列未使用
- **File**: `config/samples.tsv`
- **Problem**: `Directory` 列存在但 workflow 中 `link_fastq`/`link_blow5` 使用 `data/{sample}/` 模式，不读取此列。.test 版本甚至没有此列。
- **Fix**: 移除 Directory 列，或在 link 规则中实际使用它

### 20. `baleen_mod.py` 硬编码参数
- **File**: `workflow/scripts/baleen_mod.py` lines 23-33
- **Problem**: `proba_method`, `min_depth`, `proba_threshold`, `dedi_method` 等参数硬编码在脚本中，用户无法通过 config 调整。
- **Fix**: 通过 `snakemake.params` 从 config 传入

### 21. `mapping_list.py` 使用过时的 subprocess 模式
- **File**: `workflow/scripts/mapping_list.py`
- **Problem**: 使用 `Popen().communicate()` 循环调用外部命令，无返回值检查，出错时临时文件不清理。
- **Fix**: 改用 `subprocess.run(check=True)` 并添加 finally 块清理临时文件

### 22. `drummer.py` 硬编码工作目录
- **File**: `workflow/scripts/drummer.py` line 24
- **Problem**: `os.chdir('/opt/DRUMMER')` 假设特定容器布局。
- **Fix**: 使用相对于 snakemake 路径的动态路径

---

## Completed

- [x] 所有规则 `mem_mb` 最低 10GB — 12 个文件共 30 处已修改
- [x] `epinano_prep.sh` 诊断改进 — 捕获 stdout+stderr、扩大搜索范围、使用 params.extra
