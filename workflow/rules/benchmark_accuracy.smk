def get_all_result_tsvs(wildcards):
    """Collect all *_results.tsv paths for activated tools.

    Uses PER_COMPARISON_TOOLS and PER_SAMPLE_TOOLS from common.smk to avoid
    maintaining a separate hardcoded list.
    """
    tools = [t for t in config["tools"] if config["tools"][t]["activate"]]
    result_files = []
    for tool in tools:
        if tool in PER_COMPARISON_TOOLS:
            result_files += expand(
                f"{RESULT_ROOT}/modifications/{tool}/{{comp}}/{tool}_results.tsv",
                comp=comparisons,
            )
        elif tool in PER_SAMPLE_TOOLS:
            result_files += expand(
                f"{RESULT_ROOT}/modifications/{tool}/{{sample}}/{tool}_results.tsv",
                sample=list(samples.index),
            )
    return result_files


if config.get("benchmark", {}).get("truth_set", ""):

    rule benchmark_tool_coverage:
        input:
            results="{project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv",
            truth_set=config["benchmark"]["truth_set"],
        output:
            covered="{project}/results/benchmarks/coverage/{comparison}/{tool}_covered.tsv",
        params:
            window=config["benchmark"]["window"],
        resources:
            mem_mb=1024 * 10,
        threads: 1
        priority: 30
        log:
            "logs/{project}/benchmark_tool_coverage/{comparison}/{tool}.log",
        benchmark:
            "benchmarks/{project}/benchmark_tool_coverage/{comparison}/{tool}.benchmark.txt"
        container:
            get_container("python3")
        script:
            "../scripts/benchmark_coverage.py"


    rule benchmark_coverage_union:
        input:
            covered=lambda wc: expand(
                f"{{project}}/results/benchmarks/coverage/{{comparison}}/{{tool}}_covered.tsv",
                project=wc.project,
                comparison=wc.comparison,
                tool=get_active_comparison_tools(),
            ),
        output:
            union="{project}/results/benchmarks/coverage/{comparison}/union.tsv",
            called_sites="{project}/results/benchmarks/coverage/{comparison}/called_sites.tsv",
        resources:
            mem_mb=1024 * 10,
        threads: 1
        priority: 31
        log:
            "logs/{project}/benchmark_coverage_union/{comparison}.log",
        benchmark:
            "benchmarks/{project}/benchmark_coverage_union/{comparison}.benchmark.txt"
        container:
            get_container("python3")
        script:
            "../scripts/benchmark_coverage.py"


    rule benchmark_prediction_union:
        """Union of ALL positions reported by ANY tool.

        Unlike the coverage union (truth-covered sites only), this includes
        every position any tool reports — giving fair mode the negatives it
        needs for meaningful metrics.
        """
        input:
            results=lambda wc: expand(
                f"{{project}}/results/modifications/{{tool}}/{{comparison}}/{{tool}}_results.tsv",
                project=wc.project,
                comparison=wc.comparison,
                tool=get_active_comparison_tools(),
            ),
        output:
            prediction_union="{project}/results/benchmarks/coverage/{comparison}/union_predictions.tsv",
        resources:
            mem_mb=1024 * 10,
        threads: 1
        priority: 31
        log:
            "logs/{project}/benchmark_prediction_union/{comparison}.log",
        benchmark:
            "benchmarks/{project}/benchmark_prediction_union/{comparison}.benchmark.txt"
        container:
            get_container("python3")
        script:
            "../scripts/benchmark_coverage.py"


    rule benchmark_tool_native:
        input:
            results="{project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv",
            truth_set=config["benchmark"]["truth_set"],
        output:
            scores="{project}/results/benchmarks/per_tool/{tool}/{comparison}/native/score_comparison.tsv",
            best_metrics="{project}/results/benchmarks/per_tool/{tool}/{comparison}/native/best_metrics.tsv",
            best_score="{project}/results/benchmarks/per_tool/{tool}/{comparison}/native/best_score.txt",
        params:
            window=config["benchmark"]["window"],
        resources:
            mem_mb=1024 * 10,
        threads: 1
        priority: 32
        log:
            "logs/{project}/benchmark_tool_native/{tool}/{comparison}.log",
        benchmark:
            "benchmarks/{project}/benchmark_tool_native/{tool}/{comparison}.benchmark.txt"
        container:
            get_container("python3")
        script:
            "../scripts/benchmark_per_tool.py"


    rule benchmark_tool_fair:
        input:
            results="{project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv",
            union="{project}/results/benchmarks/coverage/{comparison}/union_predictions.tsv",
            truth_set=config["benchmark"]["truth_set"],
        output:
            scores="{project}/results/benchmarks/per_tool/{tool}/{comparison}/fair/score_comparison.tsv",
            best_metrics="{project}/results/benchmarks/per_tool/{tool}/{comparison}/fair/best_metrics.tsv",
            best_score="{project}/results/benchmarks/per_tool/{tool}/{comparison}/fair/best_score.txt",
        params:
            window=config["benchmark"]["window"],
        resources:
            mem_mb=1024 * 10,
        threads: 1
        priority: 33
        log:
            "logs/{project}/benchmark_tool_fair/{tool}/{comparison}.log",
        benchmark:
            "benchmarks/{project}/benchmark_tool_fair/{tool}/{comparison}.benchmark.txt"
        container:
            get_container("python3")
        script:
            "../scripts/benchmark_per_tool.py"


    rule benchmark_per_sample_native:
        """Benchmark per-sample tools (tandemmod, directrm, etc.) against truth set."""
        input:
            results="{project}/results/modifications/{tool}/{sample}/{tool}_results.tsv",
            truth_set=config["benchmark"]["truth_set"],
        output:
            scores="{project}/results/benchmarks/per_tool/{tool}/{sample}/native/score_comparison.tsv",
            best_metrics="{project}/results/benchmarks/per_tool/{tool}/{sample}/native/best_metrics.tsv",
            best_score="{project}/results/benchmarks/per_tool/{tool}/{sample}/native/best_score.txt",
        params:
            window=config["benchmark"]["window"],
        resources:
            mem_mb=1024 * 10,
        threads: 1
        priority: 32
        log:
            "logs/{project}/benchmark_per_sample_native/{tool}/{sample}.log",
        benchmark:
            "benchmarks/{project}/benchmark_per_sample_native/{tool}/{sample}.benchmark.txt"
        container:
            get_container("python3")
        script:
            "../scripts/benchmark_per_tool.py"


    rule benchmark_kmer_ref_negatives:
        """Generate reference-based 5-mer matched negative sites from transcriptome."""
        input:
            truth_set=config["benchmark"]["truth_set"],
            transcriptome_fasta=config["reference"]["transcriptome_fasta"],
        output:
            negatives="{project}/results/benchmarks/kmer_ref/negative_sites.tsv",
        params:
            window=max(config["benchmark"]["window"]) if isinstance(config["benchmark"]["window"], list) else config["benchmark"]["window"],
        resources:
            mem_mb=1024 * 10,
        threads: 1
        priority: 31
        log:
            "logs/{project}/benchmark_kmer_ref_negatives.log",
        benchmark:
            "benchmarks/{project}/benchmark_kmer_ref_negatives.benchmark.txt"
        conda:
            "../envs/pyfastx.yaml"
        container:
            get_container("python3")
        script:
            "../scripts/benchmark_kmer_ref_negatives.py"


    rule benchmark_tool_kmer_ref_native:
        """Per-tool benchmark using kmer-ref negatives, native eval (tool's own predictions only)."""
        input:
            results="{project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv",
            truth_set=config["benchmark"]["truth_set"],
            negatives="{project}/results/benchmarks/kmer_ref/negative_sites.tsv",
        output:
            scores="{project}/results/benchmarks/per_tool/{tool}/{comparison}/kmer_ref_native/score_comparison.tsv",
            best_metrics="{project}/results/benchmarks/per_tool/{tool}/{comparison}/kmer_ref_native/best_metrics.tsv",
            best_score="{project}/results/benchmarks/per_tool/{tool}/{comparison}/kmer_ref_native/best_score.txt",
        params:
            window=config["benchmark"]["window"],
        resources:
            mem_mb=1024 * 10,
        threads: 1
        priority: 32
        log:
            "logs/{project}/benchmark_tool_kmer_ref_native/{tool}/{comparison}.log",
        benchmark:
            "benchmarks/{project}/benchmark_tool_kmer_ref_native/{tool}/{comparison}.benchmark.txt"
        container:
            get_container("python3")
        script:
            "../scripts/benchmark_per_tool.py"


    rule benchmark_tool_kmer_ref_fair:
        """Per-tool benchmark using kmer-ref negatives, fair eval (union of all tools' predictions)."""
        input:
            results="{project}/results/modifications/{tool}/{comparison}/{tool}_results.tsv",
            union="{project}/results/benchmarks/coverage/{comparison}/union_predictions.tsv",
            truth_set=config["benchmark"]["truth_set"],
            negatives="{project}/results/benchmarks/kmer_ref/negative_sites.tsv",
        output:
            scores="{project}/results/benchmarks/per_tool/{tool}/{comparison}/kmer_ref_fair/score_comparison.tsv",
            best_metrics="{project}/results/benchmarks/per_tool/{tool}/{comparison}/kmer_ref_fair/best_metrics.tsv",
            best_score="{project}/results/benchmarks/per_tool/{tool}/{comparison}/kmer_ref_fair/best_score.txt",
        params:
            window=config["benchmark"]["window"],
        resources:
            mem_mb=1024 * 10,
        threads: 1
        priority: 33
        log:
            "logs/{project}/benchmark_tool_kmer_ref_fair/{tool}/{comparison}.log",
        benchmark:
            "benchmarks/{project}/benchmark_tool_kmer_ref_fair/{tool}/{comparison}.benchmark.txt"
        container:
            get_container("python3")
        script:
            "../scripts/benchmark_per_tool.py"


    rule benchmark_aggregate:
        input:
            native_metrics=lambda wc: (
                expand(
                    f"{{project}}/results/benchmarks/per_tool/{{tool}}/{{comparison}}/native/best_metrics.tsv",
                    project=wc.project,
                    tool=get_active_comparison_tools(),
                    comparison=comparisons,
                ) + expand(
                    f"{{project}}/results/benchmarks/per_tool/{{tool}}/{{sample}}/native/best_metrics.tsv",
                    project=wc.project,
                    tool=get_active_per_sample_tools(),
                    sample=list(samples.index),
                )
            ),
            native_all_scores=lambda wc: (
                expand(
                    f"{{project}}/results/benchmarks/per_tool/{{tool}}/{{comparison}}/native/score_comparison.tsv",
                    project=wc.project,
                    tool=get_active_comparison_tools(),
                    comparison=comparisons,
                ) + expand(
                    f"{{project}}/results/benchmarks/per_tool/{{tool}}/{{sample}}/native/score_comparison.tsv",
                    project=wc.project,
                    tool=get_active_per_sample_tools(),
                    sample=list(samples.index),
                )
            ),
            fair_metrics=lambda wc: expand(
                f"{{project}}/results/benchmarks/per_tool/{{tool}}/{{comparison}}/fair/best_metrics.tsv",
                project=wc.project,
                tool=get_active_comparison_tools(),
                comparison=comparisons,
            ),
            fair_all_scores=lambda wc: expand(
                f"{{project}}/results/benchmarks/per_tool/{{tool}}/{{comparison}}/fair/score_comparison.tsv",
                project=wc.project,
                tool=get_active_comparison_tools(),
                comparison=comparisons,
            ),
            kmer_ref_native_metrics=lambda wc: expand(
                f"{{project}}/results/benchmarks/per_tool/{{tool}}/{{comparison}}/kmer_ref_native/best_metrics.tsv",
                project=wc.project,
                tool=get_active_comparison_tools(),
                comparison=comparisons,
            ),
            kmer_ref_native_all_scores=lambda wc: expand(
                f"{{project}}/results/benchmarks/per_tool/{{tool}}/{{comparison}}/kmer_ref_native/score_comparison.tsv",
                project=wc.project,
                tool=get_active_comparison_tools(),
                comparison=comparisons,
            ),
            kmer_ref_fair_metrics=lambda wc: expand(
                f"{{project}}/results/benchmarks/per_tool/{{tool}}/{{comparison}}/kmer_ref_fair/best_metrics.tsv",
                project=wc.project,
                tool=get_active_comparison_tools(),
                comparison=comparisons,
            ),
            kmer_ref_fair_all_scores=lambda wc: expand(
                f"{{project}}/results/benchmarks/per_tool/{{tool}}/{{comparison}}/kmer_ref_fair/score_comparison.tsv",
                project=wc.project,
                tool=get_active_comparison_tools(),
                comparison=comparisons,
            ),
            called_sites=lambda wc: expand(
                f"{{project}}/results/benchmarks/coverage/{{comparison}}/called_sites.tsv",
                project=wc.project,
                comparison=comparisons,
            ),
            truth_set=config["benchmark"]["truth_set"],
        output:
            summary="{project}/results/benchmarks/cross_tool/accuracy_summary.tsv",
            overall="{project}/results/benchmarks/cross_tool/accuracy_summary_overall.tsv",
            by_comparison="{project}/results/benchmarks/cross_tool/accuracy_summary_by_comparison.tsv",
            by_negative_type="{project}/results/benchmarks/cross_tool/accuracy_summary_by_negative_type.tsv",
            by_tool="{project}/results/benchmarks/cross_tool/by_tool.tsv",
            best_scores="{project}/results/benchmarks/cross_tool/best_scores.tsv",
            called_sites_comp="{project}/results/benchmarks/cross_tool/called_sites_by_comparison.tsv",
            called_sites_sum="{project}/results/benchmarks/cross_tool/called_sites_summary.tsv",
            done=touch("{project}/results/benchmarks/.benchmark_complete"),
        resources:
            mem_mb=1024 * 10,
        threads: 1
        priority: 34
        log:
            "logs/{project}/benchmark_aggregate/aggregate.log",
        benchmark:
            "benchmarks/{project}/benchmark_aggregate/aggregate.benchmark.txt"
        container:
            get_container("python3")
        script:
            "../scripts/benchmark_aggregation.py"
