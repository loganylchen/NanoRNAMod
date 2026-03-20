rule aggregate_benchmarks:
    input:
        expand(
            "benchmarks/{project}/{stem}.benchmark.txt",
            project=PROJECT,
            stem=glob_wildcards(f"benchmarks/{PROJECT}/{{stem}}.benchmark.txt").stem,
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
    script:
        "../scripts/aggregate_benchmarks.py"
