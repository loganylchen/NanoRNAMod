# Benchmark Visualization Rules
# Generate plots and HTML reports for benchmark results

import os


rule benchmark_visualization:
    """Generate PR curves, ROC curves, and HTML benchmark report."""
    input:
        accuracy="{project}/results/benchmarks/accuracy_summary.tsv",
        resource="{project}/results/benchmarks/resource_summary.tsv",
    output:
        html="{project}/results/benchmarks/viz/benchmark_report.html",
        done=touch("{project}/results/benchmarks/viz/.done"),
    params:
        window=lambda wc: config["benchmark"].get("window", [0]),
        output_dir="{project}/results/benchmarks/viz",
    resources:
        mem_mb=1024 * 4,
    threads: 1
    priority: 25
    log:
        "logs/{project}/benchmark_visualization/viz.log",
    conda:
        "../envs/matplotlib.yaml"
    container:
        None  # Disable singularity for this rule - uses host conda
    script:
        "../scripts/benchmark_plots.py"


rule benchmark_threshold:
    """Find optimal score thresholds for each tool."""
    input:
        results=lambda wc: get_all_result_tsvs(wc),
        truth_set=config["benchmark"]["truth_set"],
    output:
        "{project}/results/benchmarks/optimal_thresholds.tsv",
    params:
        window=config["benchmark"]["window"],
        criterion="f1",
    resources:
        mem_mb=1024 * 8,
    threads: 1
    priority: 20
    log:
        "logs/{project}/benchmark_threshold/threshold.log",
    benchmark:
        "benchmarks/{project}/benchmark_threshold.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    container:
        None  # Disable singularity for this rule - uses host conda
    script:
        "../scripts/benchmark_threshold.py"


rule benchmark_detailed:
    """Generate detailed site-by-site comparison reports."""
    input:
        results=lambda wc: get_all_result_tsvs(wc),
        truth_set=config["benchmark"]["truth_set"],
    output:
        predictions="{project}/results/benchmarks/detailed_predictions.tsv",
        truth="{project}/results/benchmarks/detailed_truth.tsv",
    params:
        window=config["benchmark"]["window"],
    resources:
        mem_mb=1024 * 8,
    threads: 1
    priority: 20
    log:
        "logs/{project}/benchmark_detailed/detailed.log",
    benchmark:
        "benchmarks/{project}/benchmark_detailed.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    container:
        None  # Disable singularity for this rule - uses host conda
    script:
        "../scripts/benchmark_detailed.py"


rule benchmark_multithreshold:
    """Multi-threshold evaluation with automatic threshold generation."""
    input:
        results=lambda wc: get_all_result_tsvs(wc),
        truth_set=config["benchmark"]["truth_set"],
    output:
        evaluation="{project}/results/benchmarks/threshold_evaluation.tsv",
        optimal="{project}/results/benchmarks/optimal_thresholds_detailed.tsv",
        distribution="{project}/results/benchmarks/score_distributions.tsv",
    params:
        window=config["benchmark"]["window"],
        n_thresholds=config.get("benchmark", {}).get("n_thresholds", 50),
        custom_thresholds=config.get("benchmark", {}).get("custom_thresholds", []),
    resources:
        mem_mb=1024 * 16,
    threads: get_threads("benchmark_multithreshold", 2)
    priority: 15
    log:
        "logs/{project}/benchmark_multithreshold/multithreshold.log",
    benchmark:
        "benchmarks/{project}/benchmark_multithreshold.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    container:
        None  # Disable singularity for this rule - uses host conda
    script:
        "../scripts/benchmark_multithreshold.py"


rule benchmark_score_optimization:
    """Per-tool internal score column optimization - find best score column and threshold for each tool."""
    input:
        results=lambda wc: get_all_result_tsvs(wc),
        truth_set=config["benchmark"]["truth_set"],
    output:
        optimal="{project}/results/benchmarks/optimal_score_per_tool.tsv",
        all_eval="{project}/results/benchmarks/all_scores_evaluation.tsv",
    params:
        window=config["benchmark"]["window"],
        n_thresholds=config.get("benchmark", {}).get("n_thresholds", 100),
    resources:
        mem_mb=1024 * 16,
    threads: get_threads("benchmark_score_optimization", 2)
    priority: 16
    log:
        "logs/{project}/benchmark_score_optimization/score_optimization.log",
    benchmark:
        "benchmarks/{project}/benchmark_score_optimization.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    container:
        None  # Disable singularity for this rule - uses host conda
    script:
        "../scripts/benchmark_score_optimization.py"


rule benchmark_r_figures:
    """Generate Nature-quality figures using R ggplot2."""
    input:
        accuracy="{project}/results/benchmarks/accuracy_summary.tsv",
        accuracy_overall="{project}/results/benchmarks/accuracy_summary_overall.tsv",
        optimal_scores="{project}/results/benchmarks/optimal_score_per_tool.tsv",
        thresholds="{project}/results/benchmarks/threshold_evaluation.tsv",
        resources="{project}/results/benchmarks/resource_summary.tsv",
    output:
        dir=directory("{project}/results/benchmarks/figures"),
    params:
        window=config["benchmark"]["window"],
    resources:
        mem_mb=1024 * 8,
    threads: 1
    priority: 30
    log:
        "logs/{project}/benchmark_r_figures/figures.log",
    conda:
        "../envs/r_viz.yaml"
    container:
        None  # Disable singularity for this rule - uses host conda
    script:
        "../scripts/R/run_all_figures.R"
