rule aggregate_benchmarks:
    input:
        # Depend on tool result TSVs (not benchmark files directly — those are
        # side-effects written by Snakemake and may not exist until jobs finish).
        # By the time all result TSVs exist, all benchmark files are guaranteed
        # to exist as well, so the script's glob over benchmark_dir is safe.
        tool_results=lambda wc: get_all_result_tsvs(wc),
    output:
        "{project}/results/benchmarks/resource_summary.tsv",
    params:
        benchmark_dir="benchmarks/{project}",
    resources:
        mem_mb=1024 * 10,
    threads: 1
    priority: 30
    log:
        "logs/{project}/aggregate_benchmarks/aggregate.log",
    container:
        get_container("python3")
    script:
        "../scripts/aggregate_benchmarks.py"
