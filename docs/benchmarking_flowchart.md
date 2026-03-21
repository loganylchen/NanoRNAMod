# NanoRNAMod Benchmarking Module Flowchart

## 1. Overall Workflow Architecture

```mermaid
flowchart TB
    subgraph Input["📥 Input Data"]
        FASTQ["data/{sample}/fastq/pass.fq.gz"]
        BLOW5["data/{sample}/blow5/nanopore.blow5"]
        REF["Reference Files<br/>(ref.fa, ref.gtf)"]
        SAMPLES["samples.tsv<br/>(Native + Control)"]
        TRUTH["Truth Set<br/>(benchmark.truth_set)"]
    end

    subgraph PrepPhase["🔧 Phase 1: Data Preparation"]
        link_fastq["link_fastq<br/>Symlink FASTQ"]
        link_blow5["link_blow5<br/>Symlink BLOW5"]
        prep_ref["prep_reference<br/>Index & Process"]
        minimap2["minimap2_align<br/>Genome + Transcriptome"]
        filter["filter_reads<br/>Quality Filtering"]
        f5c_index["f5c index<br/>Signal Index"]
        f5c_event["f5c eventalign<br/>Signal Alignment"]
    end

    subgraph ModetectPhase["🔬 Phase 2: Modification Detection"]
        direction TB

        subgraph PerComparison["Per-Comparison Tools"]
            xpore["xpore"]
            nanocompore["nanocompore"]
            baleen["baleen"]
            pybaleen["pybaleen"]
            differr["differr"]
            drummer["drummer"]
            eligos2["eligos2"]
            epinano["epinano"]
            psipore["psipore"]
        end

        subgraph PerSample["Per-Sample Tools"]
            tandemmod["tandemmod"]
            directrm["directrm"]
            m6atm["m6atm"]
            rnano["rnano"]
            nanopsu["nanopsu"]
            nanomud["nanomud"]
            penguin["penguin"]
        end
    end

    subgraph PostPhase["📊 Phase 3: Post-Processing"]
        postprocess["Post-processing<br/>Format Normalization"]
        annotate["GTF Annotation<br/>(if transcriptome_gtf)"]
    end

    subgraph BenchmarkPhase["📈 Phase 4: Benchmarking"]
        resource_agg["aggregate_benchmarks<br/>Resource Summary"]
        accuracy["accuracy_benchmark<br/>Precision/Recall/F1/etc."]
    end

    subgraph Output["📤 Output"]
        results["results/modifications/{tool}/<br/>{tool}_results.tsv"]
        annotated["results/modifications/{tool}/<br/>{tool}_annotated_results.tsv"]
        resource_summary["results/benchmarks/<br/>resource_summary.tsv"]
        accuracy_summary["results/benchmarks/<br/>accuracy_summary.tsv"]
        accuracy_overall["results/benchmarks/<br/>accuracy_summary_overall.tsv"]
    end

    %% Input connections
    FASTQ --> link_fastq
    BLOW5 --> link_blow5
    REF --> prep_ref
    SAMPLES -.->|defines comparisons| ModetectPhase

    %% Prep phase connections
    link_fastq --> minimap2
    link_blow5 --> f5c_index
    prep_ref --> minimap2
    minimap2 --> filter
    filter --> f5c_index
    f5c_index --> f5c_event

    %% To detection tools
    f5c_event --> PerComparison
    f5c_event --> PerSample
    minimap2 --> PerComparison
    minimap2 --> PerSample

    %% Detection to post-processing
    PerComparison --> postprocess
    PerSample --> postprocess
    postprocess --> results
    postprocess --> annotate
    annotate --> annotated

    %% To benchmarking
    results -.->|all *_results.tsv| accuracy
    TRUTH -.->|ground truth| accuracy
    results --> resource_agg

    %% Benchmarking outputs
    resource_agg --> resource_summary
    accuracy --> accuracy_summary
    accuracy --> accuracy_overall

    %% Styling
    classDef inputStyle fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    classDef prepStyle fill:#fff3e0,stroke:#e65100,stroke-width:2px
    classDef detectStyle fill:#f3e5f5,stroke:#4a148c,stroke-width:2px
    classDef benchStyle fill:#e8f5e9,stroke:#1b5e20,stroke-width:2px
    classDef outputStyle fill:#fce4ec,stroke:#880e4f,stroke-width:2px

    class Input inputStyle
    class PrepPhase prepStyle
    class ModetectPhase detectStyle
    class BenchmarkPhase benchStyle
    class Output outputStyle
```

## 2. Benchmarking Module Detail

```mermaid
flowchart LR
    subgraph Inputs
        tool_results["Tool Results<br/>*_results.tsv"]
        truth["Truth Set<br/>(transcript, position,<br/>modification_type, label)"]
    end

    subgraph Processing["accuracy_benchmark.py"]
        normalize["Column Normalization<br/>(transcript, position)"]
        detect_score["Score Column Detection<br/>(tool-specific)"]
        match["Position Matching<br/>(within window)"]
        metrics["Metrics Computation"]
    end

    subgraph Metrics["Output Metrics"]
        m1["Precision"]
        m2["Recall"]
        m3["F1 Score"]
        m4["Specificity"]
        m5["MCC"]
        m6["AUPRC"]
        m7["AUROC"]
        m8["Confusion Matrix<br/>(TP/FP/FN/TN)"]
    end

    tool_results --> normalize
    truth --> match
    normalize --> detect_score
    detect_score --> match
    match --> metrics
    metrics --> m1 & m2 & m3 & m4 & m5 & m6 & m7 & m8
```

## 3. Tool Categories & Data Flow

```mermaid
flowchart TB
    subgraph ComparisonTools["Per-Comparison Tools<br/>(Native vs Control)"]
        C1["xpore"]
        C2["nanocompore"]
        C3["baleen / pybaleen"]
        C4["differr"]
        C5["drummer"]
        C6["eligos2"]
        C7["epinano"]
        C8["psipore"]
    end

    subgraph SampleTools["Per-Sample Tools<br/>(Single Sample)"]
        S1["tandemmod"]
        S2["directrm"]
        S3["m6atm"]
        S4["rnano"]
        S5["nanopsu"]
        S6["nanomud"]
        S7["penguin"]
    end

    subgraph Outputs
        O1["results/modifications/{tool}/{comp}/<br/>{tool}_results.tsv"]
        O2["results/modifications/{tool}/{sample}/<br/>{tool}_results.tsv"]
    end

    ComparisonTools --> O1
    SampleTools --> O2

    O1 & O2 --> B["accuracy_benchmark"]
```

## 4. Resource Benchmarking Flow

```mermaid
flowchart LR
    subgraph Snakemake["Snakemake Built-in Benchmarking"]
        rule["Each Rule Execution"]
        bench_file["benchmarks/{project}/{stem}.benchmark.txt"]
    end

    subgraph Metrics["Resource Metrics"]
        time["Wall-clock Time<br/>(s, h:m:s)"]
        memory["Memory Usage<br/>(max_rss, max_vms, max_uss, max_pss)"]
        io["I/O Operations<br/>(io_in, io_out)"]
        cpu["CPU Metrics<br/>(mean_load, cpu_time)"]
    end

    subgraph Aggregate["aggregate_benchmarks.py"]
        collect["Collect All .benchmark.txt"]
        summarize["Summarize by Tool & Sample"]
    end

    rule --> bench_file
    bench_file --> time & memory & io & cpu
    time & memory & io & cpu --> collect
    collect --> summarize
    summarize --> output["resource_summary.tsv"]
```

---

## 5. Visualization Design Recommendations

### 5.1 Recommended Tools

| Tool | Purpose | Reason |
|------|---------|--------|
| **Matplotlib + Seaborn** | Static charts | Python-native, integrates with workflow |
| **Plotly** | Interactive charts | Hover tooltips, zoom, export |
| **MultiQC** | QC aggregation | Already in bioinformatics ecosystem |
| **Vega-Lite / Altair** | Declarative viz | JSON export, Jupyter integration |

### 5.2 Color Palette

```
Primary Tools:
- xpore:        #1f77b4 (blue)
- nanocompore:  #ff7f0e (orange)
- baleen:       #2ca02c (green)
- pybaleen:     #17becf (cyan)
- differr:      #d62728 (red)
- drummer:      #9467bd (purple)
- eligos2:      #8c564b (brown)
- epinano:      #e377c2 (pink)

Metric Categories:
- Accuracy:     #2196f3 (blue)
- Performance:  #4caf50 (green)
- Resource:     #ff9800 (orange)
- Quality:      #9c27b0 (purple)
```

### 5.3 Suggested Visualizations

#### A. Accuracy Comparison Dashboard

```
┌─────────────────────────────────────────────────────────────┐
│  F1 Score by Tool (Bar Chart)                               │
│  ┌──────────────────────────────────────────────────────┐   │
│  │ ████████████████████████████  xpore (0.85)           │   │
│  │ ██████████████████████  nanocompore (0.78)           │   │
│  │ ███████████████████  baleen (0.72)                   │   │
│  └──────────────────────────────────────────────────────┘   │
├─────────────────────────────────────────────────────────────┤
│  Precision vs Recall (Scatter Plot)                         │
│  ┌──────────────────────────────────────────────────────┐   │
│  │      R                                                │   │
│  │      e    ● xpore                                     │   │
│  │      c       ● nano                                   │   │
│  │      a          ● baleen                              │   │
│  │      l    ● pybaleen                                  │   │
│  │             Precision →                               │   │
│  └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

#### B. Resource Usage Heatmap

```
                    Runtime (h)    Memory (GB)    CPU (%)
    xpore           [███░░░] 2.1   [██████] 8.2   [████░] 75%
    nanocompore     [████░░] 3.5   [█████░] 6.5   [█████] 85%
    baleen          [██░░░░] 1.2   [██░░░░] 3.1   [██░░░] 45%
    pybaleen        [██░░░░] 1.0   [█░░░░░] 2.0   [███░░] 55%
```

#### C. ROC/PR Curves (Per Tool)

```python
# Pseudo-code for generating ROC curves
for tool in tools:
    fpr, tpr, thresholds = roc_curve(y_true, y_scores[tool])
    plt.plot(fpr, tpr, label=f'{tool} (AUC = {auroc[tool]:.2f})')
```

### 5.4 Output File Structure

```
results/benchmarks/
├── accuracy_summary.tsv           # Per modification_type metrics
├── accuracy_summary_overall.tsv   # Aggregated metrics
├── resource_summary.tsv           # Runtime & memory usage
├── plots/
│   ├── f1_comparison.png          # Bar chart
│   ├── precision_recall.png       # Scatter plot
│   ├── roc_curves.png             # Multi-line plot
│   ├── resource_heatmap.png       # Resource comparison
│   └── per_tool/
│       ├── xpore_roc.png
│       ├── xpore_pr.png
│       └── ...
└── interactive/
    └── benchmark_dashboard.html   # Plotly dashboard
```

---

## 6. Implementation Notes

### 6.1 Key Dependencies

- `sklearn.metrics` for AUROC/AUPRC computation
- `pandas` for data manipulation
- `matplotlib`/`seaborn` for static plots
- `plotly` for interactive dashboards (optional)

### 6.2 Configuration

```yaml
# config/config.yaml
benchmark:
  truth_set: "path/to/truth_set.tsv"
  window: [0, 5, 10]  # Multi-window evaluation

  visualization:
    generate_plots: true
    interactive: true
    output_format: ["png", "html"]
```

### 6.3 Future Enhancements

1. **Bootstrap confidence intervals** for metric uncertainty
2. **Stratified analysis** by transcript region (5'UTR, CDS, 3'UTR)
3. **Per-modification-type breakdown** (m6A vs Psi vs m1A)
4. **Sample size sensitivity** analysis
5. **Tool consensus** visualization (Venn diagrams, upset plots)
