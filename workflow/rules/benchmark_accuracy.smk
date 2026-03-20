def get_all_result_tsvs(wildcards):
    """Collect all *_results.tsv paths for activated tools."""
    tools = [t for t in config["tools"] if config["tools"][t]["activate"]]
    result_files = []
    per_comparison_tools = ["xpore", "nanocompore", "baleen", "differr", "drummer",
                             "eligos2", "epinano", "psipore"]
    per_sample_tools = ["tandemmod", "directrm", "m6atm", "rnano",
                        "nanopsu", "nanomud", "penguin"]
    for tool in tools:
        if tool in per_comparison_tools:
            result_files += expand(
                f"{RESULT_ROOT}/modifications/{tool}/{{comp}}/{tool}_results.tsv",
                comp=comparisons,
            )
        elif tool in per_sample_tools:
            result_files += expand(
                f"{RESULT_ROOT}/modifications/{tool}/{{sample}}/{tool}_results.tsv",
                sample=list(samples.index),
            )
    return result_files


if config.get("benchmark", {}).get("truth_set", ""):
    rule accuracy_benchmark:
        input:
            results=get_all_result_tsvs,
            truth_set=config["benchmark"]["truth_set"],
        output:
            expand("{project}/results/benchmarks/accuracy_summary{suffix}.tsv", suffix=["", "_overall"]),
        params:
            window=config["benchmark"]["window"],
        resources:
            mem_mb=1024 * 8,
        threads: 1
        priority: 30
        log:
            "logs/{project}/accuracy_benchmark/accuracy.log",
        conda:
            "../envs/pandas.yaml"
        script:
            "../scripts/accuracy_benchmark.py"
