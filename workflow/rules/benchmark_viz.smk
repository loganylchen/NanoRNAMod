# Benchmark Visualization Rules
# Generate plots and HTML reports for benchmark results

import os


rule benchmark_visualization:
    """Generate PR curves, ROC curves, and HTML benchmark report."""
    input:
        # Depend on the touch file from accuracy_benchmark
        # This ensures the accuracy_benchmark rule has completed
        done="{project}/results/benchmarks/.benchmark_complete",
        accuracy="{project}/results/benchmarks/cross_tool/accuracy_summary.tsv",
        resource="{project}/results/benchmarks/resource_summary.tsv",
    output:
        html="{project}/results/benchmarks/viz/benchmark_report.html",
        done_viz=touch("{project}/results/benchmarks/viz/.done"),
    params:
        window=lambda wc: config["benchmark"].get("window", [0]),
        output_dir="{project}/results/benchmarks/viz",
    resources:
        mem_mb=1024 * 10,
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
        mem_mb=1024 * 10,
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
        mem_mb=1024 * 10,
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
        mem_mb=1024 * 10,
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


rule benchmark_same_base_negatives:
    """Same-base negative benchmarking - use sites with same nucleotide as negative controls (Strategy 2)."""
    input:
        results=lambda wc: get_all_result_tsvs(wc),
        truth_set=config["benchmark"]["truth_set"],
    output:
        metrics="{project}/results/benchmarks/accuracy_summary_same_base_negatives.tsv",
        negatives="{project}/results/benchmarks/same_base_negative_sites.tsv",
    params:
        window=config["benchmark"]["window"],
        base_column=config.get("benchmark", {}).get("base_column", "auto"),
    resources:
        mem_mb=1024 * 10,
    threads: 1
    priority: 21
    log:
        "logs/{project}/benchmark_same_base_negatives/same_base_negatives.log",
    benchmark:
        "benchmarks/{project}/benchmark_same_base_negatives.benchmark.txt"
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_same_base_negatives.py"


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
        score_summary="{project}/results/benchmarks/score_column_summary.tsv",
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


rule benchmark_statistics:
    """Compute bootstrap confidence intervals, significance tests, and effect sizes.

    Runs AFTER accuracy_benchmark to ensure metrics are available.
    Computes 95% bootstrap CIs for all metrics, paired Wilcoxon tests,
    permutation tests, FDR correction, and effect sizes (Cohen's d, Cliff's delta).
    """
    input:
        summary="{project}/results/benchmarks/cross_tool/accuracy_summary.tsv",
        by_comparison="{project}/results/benchmarks/cross_tool/accuracy_summary_by_comparison.tsv",
    output:
        ci="{project}/results/benchmarks/statistics/bootstrap_ci.tsv",
        significance="{project}/results/benchmarks/statistics/significance_tests.tsv",
        fdr="{project}/results/benchmarks/statistics/fdr_corrected.tsv",
        effects="{project}/results/benchmarks/statistics/effect_sizes.tsv",
        auc_comparison="{project}/results/benchmarks/statistics/auc_comparison.tsv",
    params:
        n_bootstrap=config.get("benchmark", {}).get("n_bootstrap", 1000),
        alpha=config.get("benchmark", {}).get("alpha", 0.05),
        fdr_method=config.get("benchmark", {}).get("fdr_method", "benjamini-hochberg"),
    resources:
        mem_mb=1024 * 10,
    threads: 4
    priority: 32
    log:
        "logs/{project}/benchmark_statistics/statistics.log",
    benchmark:
        "benchmarks/{project}/benchmark_statistics.benchmark.txt"
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_statistics.py"


rule benchmark_sensitivity:
    """Sensitivity analysis for coverage depth and threshold robustness.

    Coverage data source: Uses depth information from alignment BAM files.
    The detailed_predictions.tsv must include a 'coverage' column populated
    by accuracy_benchmark.py from per-site depth calculations.

    If coverage column is not present, the analysis falls back to using
    read count as a proxy for coverage depth.
    """
    input:
        predictions="{project}/results/benchmarks/detailed_predictions.tsv",
        truth_set=config["benchmark"]["truth_set"],
    output:
        coverage="{project}/results/benchmarks/sensitivity/coverage_analysis.tsv",
        score_dist="{project}/results/benchmarks/sensitivity/score_distribution.tsv",
        threshold_robust="{project}/results/benchmarks/sensitivity/threshold_robustness.tsv",
    params:
        coverage_bins=config.get("benchmark", {}).get("coverage_bins", [0, 10, 20, 50, 100, 200, 500]),
        n_splits=config.get("benchmark", {}).get("n_splits", 5),
    resources:
        mem_mb=1024 * 10,
    threads: 2
    priority: 33
    log:
        "logs/{project}/benchmark_sensitivity/sensitivity.log",
    benchmark:
        "benchmarks/{project}/benchmark_sensitivity.benchmark.txt"
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_sensitivity.py"


rule benchmark_r_figures:
    """Generate Nature-quality figures using R ggplot2 with statistical overlays."""
    input:
        # Cross-tool accuracy metrics
        aggregated=expand("{{project}}/results/benchmarks/cross_tool/{agg_file}", agg_file=[
            "accuracy_summary.tsv",
            "accuracy_summary_overall.tsv",
            "accuracy_summary_by_comparison.tsv",
            "accuracy_summary_by_negative_type.tsv",
        ]),
        optimal_scores="{project}/results/benchmarks/cross_tool/best_scores.tsv",
        by_tool="{project}/results/benchmarks/cross_tool/by_tool.tsv",
        score_optimization="{project}/results/benchmarks/optimal_score_per_tool.tsv",
        score_column_summary="{project}/results/benchmarks/score_column_summary.tsv",
        thresholds="{project}/results/benchmarks/threshold_evaluation.tsv",
        resources="{project}/results/benchmarks/resource_summary.tsv",
        resource_by_tool="{project}/results/benchmarks/resource_by_tool.tsv",
        # Statistical analysis results
        bootstrap_ci="{project}/results/benchmarks/statistics/bootstrap_ci.tsv",
        significance="{project}/results/benchmarks/statistics/significance_tests.tsv",
        fdr="{project}/results/benchmarks/statistics/fdr_corrected.tsv",
        effect_sizes="{project}/results/benchmarks/statistics/effect_sizes.tsv",
        auc_comparison="{project}/results/benchmarks/statistics/auc_comparison.tsv",
        # Sensitivity analysis results
        coverage="{project}/results/benchmarks/sensitivity/coverage_analysis.tsv",
        score_dist="{project}/results/benchmarks/sensitivity/score_distribution.tsv",
        threshold_robust="{project}/results/benchmarks/sensitivity/threshold_robustness.tsv",
    output:
        # Main figures (flat layout with co-located source data)
        fig1="{project}/results/benchmarks/figures/fig1_overall_accuracy.pdf",
        fig2="{project}/results/benchmarks/figures/fig2_precision_recall.pdf",
        fig3="{project}/results/benchmarks/figures/fig3_auroc.pdf",
        fig4="{project}/results/benchmarks/figures/fig4_native_vs_fair.pdf",
        fig5="{project}/results/benchmarks/figures/fig5_best_score.pdf",
        fig6="{project}/results/benchmarks/figures/fig6_coverage_sensitivity.pdf",
        fig7="{project}/results/benchmarks/figures/fig7_resource_usage.pdf",
        fig8="{project}/results/benchmarks/figures/fig8_tool_ranking.pdf",
        # Supplementary figures
        sfig1="{project}/results/benchmarks/figures/sfig1_per_comparison.pdf",
        sfig2="{project}/results/benchmarks/figures/sfig2_native_vs_fair_detail.pdf",
        sfig3="{project}/results/benchmarks/figures/sfig3_threshold_robustness.pdf",
        sfig4="{project}/results/benchmarks/figures/sfig4_effect_sizes.pdf",
        sfig5="{project}/results/benchmarks/figures/sfig5_score_heatmap.pdf",
        # Source data files (co-located in figures/)
        data1="{project}/results/benchmarks/figures/fig1_overall_accuracy_data.tsv",
        data2="{project}/results/benchmarks/figures/fig2_precision_recall_data.tsv",
        data3="{project}/results/benchmarks/figures/fig3_auroc_data.tsv",
        data4="{project}/results/benchmarks/figures/fig4_native_vs_fair_data.tsv",
        data5="{project}/results/benchmarks/figures/fig5_best_score_data.tsv",
        data6="{project}/results/benchmarks/figures/fig6_coverage_sensitivity_data.tsv",
        data7="{project}/results/benchmarks/figures/fig7_resource_usage_data.tsv",
        data8="{project}/results/benchmarks/figures/fig8_tool_ranking_data.tsv",
    params:
        window=config["benchmark"]["window"],
        theme=config.get("benchmark", {}).get("figure_theme", "nature"),
        dpi=config.get("benchmark", {}).get("dpi", 300),
        format=config.get("benchmark", {}).get("fig_format", "pdf"),
    resources:
        mem_mb=1024 * 10,
    threads: 1
    priority: 34
    log:
        "logs/{project}/benchmark_r_figures/figures.log",
    container:
        get_container("r_viz")
    script:
        "../scripts/R/run_all_figures.R"


rule benchmark_pdf_report:
    """Generate comprehensive PDF report with overall, per-comparison, and per-tool analysis."""
    input:
        # Depend on touch file
        done="{project}/results/benchmarks/.benchmark_complete",
        thresholds="{project}/results/benchmarks/threshold_evaluation.tsv",
        optimal="{project}/results/benchmarks/optimal_thresholds_detailed.tsv",
        distributions="{project}/results/benchmarks/score_distributions.tsv",
    output:
        pdf="{project}/results/benchmarks/benchmark_report.pdf",
    params:
        benchmark_dir="{project}/results/benchmarks/cross_tool",
    resources:
        mem_mb=1024 * 10,
    threads: 1
    priority: 35
    log:
        "logs/{project}/benchmark_pdf_report/report.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_pdf_report.py"


rule benchmark_per_tool_comparison_figures:
    """Generate per-tool × comparison figures (ROC, PR, score dist, native vs fair)."""
    input:
        native_scores="{project}/results/benchmarks/per_tool/{tool}/{comparison}/native/score_comparison.tsv",
        fair_scores="{project}/results/benchmarks/per_tool/{tool}/{comparison}/fair/score_comparison.tsv",
        results="{project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv",
        union_predictions="{project}/results/benchmarks/coverage/{comparison}/union_predictions.tsv",
        truth_set=config["benchmark"]["truth_set"],
    output:
        roc_pdf="{project}/results/benchmarks/per_tool/{tool}/{comparison}/figures/roc_curve.pdf",
        roc_data="{project}/results/benchmarks/per_tool/{tool}/{comparison}/figures/roc_curve_data.tsv",
        pr_pdf="{project}/results/benchmarks/per_tool/{tool}/{comparison}/figures/pr_curve.pdf",
        pr_data="{project}/results/benchmarks/per_tool/{tool}/{comparison}/figures/pr_curve_data.tsv",
        dist_pdf="{project}/results/benchmarks/per_tool/{tool}/{comparison}/figures/score_distribution.pdf",
        dist_data="{project}/results/benchmarks/per_tool/{tool}/{comparison}/figures/score_distribution_data.tsv",
        nvf_pdf="{project}/results/benchmarks/per_tool/{tool}/{comparison}/figures/native_vs_fair.pdf",
        nvf_data="{project}/results/benchmarks/per_tool/{tool}/{comparison}/figures/native_vs_fair_data.tsv",
    params:
        window=config["benchmark"]["window"],
    resources:
        mem_mb=1024 * 10,
    threads: 1
    priority: 35
    log:
        "logs/{project}/benchmark_per_tool_comparison_figures/{tool}/{comparison}.log",
    benchmark:
        "benchmarks/{project}/benchmark_per_tool_comparison_figures/{tool}/{comparison}.benchmark.txt"
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_per_tool_figures.py"


rule benchmark_per_tool_summary_figures:
    """Generate per-tool summary figures across all comparisons."""
    input:
        native_scores=lambda wc: expand(
            f"{{project}}/results/benchmarks/per_tool/{{tool}}/{{comparison}}/native/score_comparison.tsv",
            project=wc.project,
            tool=wc.tool,
            comparison=comparisons,
        ),
        fair_scores=lambda wc: expand(
            f"{{project}}/results/benchmarks/per_tool/{{tool}}/{{comparison}}/fair/score_comparison.tsv",
            project=wc.project,
            tool=wc.tool,
            comparison=comparisons,
        ),
    output:
        acc_pdf="{project}/results/benchmarks/per_tool/{tool}/figures/accuracy_by_comparison.pdf",
        acc_data="{project}/results/benchmarks/per_tool/{tool}/figures/accuracy_by_comparison_data.tsv",
        nvf_pdf="{project}/results/benchmarks/per_tool/{tool}/figures/native_vs_fair_summary.pdf",
        nvf_data="{project}/results/benchmarks/per_tool/{tool}/figures/native_vs_fair_summary_data.tsv",
        sc_pdf="{project}/results/benchmarks/per_tool/{tool}/figures/score_columns_comparison.pdf",
        sc_data="{project}/results/benchmarks/per_tool/{tool}/figures/score_columns_comparison_data.tsv",
        thresh_pdf="{project}/results/benchmarks/per_tool/{tool}/figures/threshold_curve.pdf",
        thresh_data="{project}/results/benchmarks/per_tool/{tool}/figures/threshold_curve_data.tsv",
    resources:
        mem_mb=1024 * 10,
    threads: 1
    priority: 35
    log:
        "logs/{project}/benchmark_per_tool_summary_figures/{tool}.log",
    benchmark:
        "benchmarks/{project}/benchmark_per_tool_summary_figures/{tool}.benchmark.txt"
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_per_tool_summary.py"


rule benchmark_per_sample_figures:
    """Generate per-sample tool figures (native only, no fair mode)."""
    input:
        native_scores="{project}/results/benchmarks/per_tool/{tool}/{sample}/native/score_comparison.tsv",
        results="{project}/results/modifications/{tool}/{sample}/{tool}_results.tsv",
        truth_set=config["benchmark"]["truth_set"],
    output:
        dist_pdf="{project}/results/benchmarks/per_tool/{tool}/{sample}/figures/score_distribution.pdf",
        dist_data="{project}/results/benchmarks/per_tool/{tool}/{sample}/figures/score_distribution_data.tsv",
    params:
        window=config["benchmark"]["window"],
    resources:
        mem_mb=1024 * 10,
    threads: 1
    priority: 35
    log:
        "logs/{project}/benchmark_per_sample_figures/{tool}/{sample}.log",
    benchmark:
        "benchmarks/{project}/benchmark_per_sample_figures/{tool}/{sample}.benchmark.txt"
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_per_sample_figures.py"


rule benchmark_cross_tool_curves:
    """Generate cross-tool overlaid ROC/PR curves, AUROC/AUPRC bar chart, and score distribution plots."""
    input:
        results=lambda wc: get_all_result_tsvs(wc),
        truth_set=config["benchmark"]["truth_set"],
        best_scores="{project}/results/benchmarks/cross_tool/best_scores.tsv",
        union_predictions=lambda wc: expand(
            f"{{project}}/results/benchmarks/coverage/{{comparison}}/union_predictions.tsv",
            project=wc.project,
            comparison=comparisons,
        ),
    output:
        roc_pdf="{project}/results/benchmarks/cross_tool/figures/cross_tool_roc.pdf",
        roc_data="{project}/results/benchmarks/cross_tool/figures/cross_tool_roc_data.tsv",
        pr_pdf="{project}/results/benchmarks/cross_tool/figures/cross_tool_pr.pdf",
        pr_data="{project}/results/benchmarks/cross_tool/figures/cross_tool_pr_data.tsv",
        bar_pdf="{project}/results/benchmarks/cross_tool/figures/cross_tool_auroc_auprc.pdf",
        bar_data="{project}/results/benchmarks/cross_tool/figures/cross_tool_auroc_auprc_data.tsv",
        dist_pdf="{project}/results/benchmarks/cross_tool/figures/cross_tool_score_distribution.pdf",
        dist_data="{project}/results/benchmarks/cross_tool/figures/cross_tool_score_distribution_data.tsv",
    params:
        window=config["benchmark"]["window"],
    resources:
        mem_mb=1024 * 10,
    threads: 1
    priority: 35
    log:
        "logs/{project}/benchmark_cross_tool_curves/curves.log",
    benchmark:
        "benchmarks/{project}/benchmark_cross_tool_curves.benchmark.txt"
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_cross_tool_curves.py"


rule benchmark_resource_by_tool:
    """Analyze resource usage by modification tool including prerequisite steps."""
    input:
        done="{project}/results/benchmarks/.benchmark_complete",
    output:
        tsv="{project}/results/benchmarks/resource_by_tool.tsv",
        pdf="{project}/results/benchmarks/resource_by_tool.pdf",
    params:
        benchmark_dir="benchmarks/{project}",
    resources:
        mem_mb=1024 * 10,
    threads: 1
    priority: 36
    log:
        "logs/{project}/benchmark_resource_by_tool/analysis.log",
    container:
        get_container("python3")
    script:
        "../scripts/benchmark_resource_by_tool.py"
