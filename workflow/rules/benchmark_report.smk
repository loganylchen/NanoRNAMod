rule aggregate_benchmarks:
    input:
        # Depend only on benchmark files from activated tools.
        # The script globs benchmark_dir at runtime so it picks up all step files;
        # these inputs just ensure the tool jobs have completed before aggregating.
        native_benchmarks=lambda wc: expand(
            "benchmarks/{project}/{comp}.{tool}.benchmark.txt",
            project=wc.project,
            comp=comparisons,
            tool=get_active_comparison_tools(),
        ),
        sample_benchmarks=lambda wc: expand(
            "benchmarks/{project}/{sample}.{tool}.benchmark.txt",
            project=wc.project,
            sample=list(samples.index),
            tool=get_active_per_sample_tools(),
        ),
    output:
        "{project}/results/benchmarks/resource_summary.tsv",
    params:
        benchmark_dir="benchmarks/{project}",
    resources:
        mem_mb=1024 * 4,
    threads: 1
    priority: 30
    log:
        "logs/{project}/aggregate_benchmarks/aggregate.log",
    container:
        get_container("python3")
    script:
        "../scripts/aggregate_benchmarks.py"
