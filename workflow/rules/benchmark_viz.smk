# Benchmark Visualization Rules
# Generate plots and HTML reports for benchmark results

import os


def get_per_tool_benchmark_files(wildcards, file_type):
    """Get per-tool benchmark files for aggregation.

    Args:
        file_type: One of 'accuracy', 'overall', 'by_comparison', 'by_negative_type'
    """
    tools = [t for t in config["tools"] if config["tools"][t]["activate"]]
    filename_map = {
        'accuracy': 'accuracy_summary.tsv',
        'overall': 'accuracy_summary_overall.tsv',
        'by_comparison': 'accuracy_summary_by_comparison.tsv',
        'by_negative_type': 'accuracy_summary_by_negative_type.tsv',
    }
    filename = filename_map[file_type]
    return [f"{wildcards.project}/results/benchmarks/{tool}/{filename}" for tool in tools]


rule aggregate_benchmark_results:
    """Aggregate per-tool benchmark results into flat files for visualization."""
    input:
        benchmark_complete="{project}/results/benchmarks/.benchmark_complete",
        accuracy=lambda wc: get_per_tool_benchmark_files(wc, 'accuracy'),
        overall=lambda wc: get_per_tool_benchmark_files(wc, 'overall'),
        by_comparison=lambda wc: get_per_tool_benchmark_files(wc, 'by_comparison'),
        by_negative_type=lambda wc: get_per_tool_benchmark_files(wc, 'by_negative_type'),
    output:
        accuracy="{project}/results/benchmarks/accuracy_summary.tsv",
        overall="{project}/results/benchmarks/accuracy_summary_overall.tsv",
        by_comparison="{project}/results/benchmarks/accuracy_summary_by_comparison.tsv",
        by_negative_type="{project}/results/benchmarks/accuracy_summary_by_negative_type.tsv",
    run:
        import pandas as pd
        import os

        # Aggregate accuracy_summary.tsv (per modification type)
        dfs = []
        for f in input.accuracy:
            if os.path.exists(f):
                df = pd.read_csv(f, sep='\t')
                dfs.append(df)
        if dfs:
            combined = pd.concat(dfs, ignore_index=True).sort_values(
                ["tool", "modification_type", "window", "f1"],
                ascending=[True, True, True, False]
            )
        else:
            combined = pd.DataFrame()
        combined.to_csv(output.accuracy, sep='\t', index=False)

        # Aggregate accuracy_summary_overall.tsv
        dfs = []
        for f in input.overall:
            if os.path.exists(f):
                df = pd.read_csv(f, sep='\t')
                dfs.append(df)
        if dfs:
            combined = pd.concat(dfs, ignore_index=True).sort_values(
                ["tool", "window", "f1"],
                ascending=[True, True, False]
            )
        else:
            combined = pd.DataFrame()
        combined.to_csv(output.overall, sep='\t', index=False)

        # Aggregate accuracy_summary_by_comparison.tsv
        dfs = []
        for f in input.by_comparison:
            if os.path.exists(f):
                df = pd.read_csv(f, sep='\t')
                dfs.append(df)
        if dfs:
            combined = pd.concat(dfs, ignore_index=True).sort_values(
                ["tool", "comparison", "window", "f1"],
                ascending=[True, True, True, False]
            )
        else:
            combined = pd.DataFrame()
        combined.to_csv(output.by_comparison, sep='\t', index=False)

        # Aggregate accuracy_summary_by_negative_type.tsv
        dfs = []
        for f in input.by_negative_type:
            if os.path.exists(f):
                df = pd.read_csv(f, sep='\t')
                dfs.append(df)
        if dfs:
            combined = pd.concat(dfs, ignore_index=True).sort_values(
                ["tool", "modification_type", "negative_type", "window", "f1"],
                ascending=[True, True, True, True, False]
            )
        else:
            combined = pd.DataFrame()
        combined.to_csv(output.by_negative_type, sep='\t', index=False)

        print(f"Aggregated benchmark results from {len(input.accuracy)} tools")


rule benchmark_visualization:
    """Generate PR curves, ROC curves, and HTML benchmark report."""
    input:
        aggregated=expand("{{project}}/results/benchmarks/{agg_file}", agg_file=[
            "accuracy_summary.tsv",
            "accuracy_summary_overall.tsv",
            "accuracy_summary_by_comparison.tsv",
            "accuracy_summary_by_negative_type.tsv",
        ]),
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
    container:
        get_container("python3")
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
    container:
        get_container("python3")
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
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_detailed.py"


rule benchmark_kmer_negatives:
    """K-mer negative benchmarking - use sites with same k-mer as negative controls."""
    input:
        results=lambda wc: get_all_result_tsvs(wc),
        truth_set=config["benchmark"]["truth_set"],
    output:
        metrics="{project}/results/benchmarks/accuracy_summary_kmer_negatives.tsv",
        negatives="{project}/results/benchmarks/kmer_negative_sites.tsv",
    params:
        window=config["benchmark"]["window"],
        kmer_column=config.get("benchmark", {}).get("kmer_column", "auto"),
    resources:
        mem_mb=1024 * 8,
    threads: 1
    priority: 22
    log:
        "logs/{project}/benchmark_kmer_negatives/kmer_negatives.log",
    benchmark:
        "benchmarks/{project}/benchmark_kmer_negatives.benchmark.txt"
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_kmer_negatives.py"


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
    container:
        get_container("python3")
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
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_score_optimization.py"


rule benchmark_r_figures:
    """Generate Nature-quality figures using R ggplot2."""
    input:
        aggregated=expand("{{project}}/results/benchmarks/{agg_file}", agg_file=[
            "accuracy_summary.tsv",
            "accuracy_summary_overall.tsv",
            "accuracy_summary_by_comparison.tsv",
            "accuracy_summary_by_negative_type.tsv",
        ]),
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
    container:
        get_container("r_viz")
    script:
        "../scripts/R/run_all_figures.R"


rule benchmark_pdf_report:
    """Generate comprehensive PDF report with overall, per-comparison, and per-tool analysis."""
    input:
        aggregated=expand("{{project}}/results/benchmarks/{agg_file}", agg_file=[
            "accuracy_summary.tsv",
            "accuracy_summary_overall.tsv",
            "accuracy_summary_by_comparison.tsv",
            "accuracy_summary_by_negative_type.tsv",
        ]),
        thresholds="{project}/results/benchmarks/threshold_evaluation.tsv",
        optimal="{project}/results/benchmarks/optimal_thresholds_detailed.tsv",
        distributions="{project}/results/benchmarks/score_distributions.tsv",
    output:
        pdf="{project}/results/benchmarks/benchmark_report.pdf",
    params:
        benchmark_dir="{project}/results/benchmarks",
    resources:
        mem_mb=1024 * 8,
    threads: 1
    priority: 35
    log:
        "logs/{project}/benchmark_pdf_report/report.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_pdf_report.py"


rule benchmark_resource_by_tool:
    """Analyze resource usage by modification tool including prerequisite steps."""
    input:
        benchmark_dir="benchmarks/{project}",
    output:
        tsv="{project}/results/benchmarks/resource_by_tool.tsv",
        pdf="{project}/results/benchmarks/resource_by_tool.pdf",
    params:
        benchmark_dir="benchmarks/{project}",
    resources:
        mem_mb=1024 * 4,
    threads: 1
    priority: 36
    log:
        "logs/{project}/benchmark_resource_by_tool/analysis.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_resource_by_tool.py"
