rule post_xpore_sampled:
    input:
        "results/xpore/{native}_{control}-{sample_size}-{n}/diffmod.table",
    output:
        "results/xpore/{native}_{control}-{sample_size}-{n}/xpore_results.tsv"
    params:
        tool="xpore",
    log:
        "logs/post_xpore_sampled_format/{native}_{control}_{sample_size}_{n}.log"
    benchmark:
        "benchmarks/{native}_{control}_{sample_size}_{n}.post_xpore_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule post_xpore:
    input:
        "results/xpore/{native}_{control}/diffmod.table",
    output:
        "results/xpore/{native}_{control}/xpore_results.tsv"
    params:
        tool="xpore",
    log:
        "logs/post_xpore_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_xpore_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule post_baleen_sampled:
    input:
        'results/baleen/{native}_{control}-{sample_size}-{n}/transcripts.csv'
    output:
        "results/baleen/{native}_{control}-{sample_size}-{n}/baleen_results.tsv"
    params:
        tool="baleen",
    log:
        "logs/post_baleen_sampled_format/{native}_{control}_{sample_size}_{n}.log"
    benchmark:
        "benchmarks/{native}_{control}_{sample_size}_{n}.post_baleen_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule post_baleen:
    input:
        'results/baleen/{native}_{control}/transcripts.csv'
    output:
        "results/baleen/{native}_{control}/baleen_results.tsv"
    params:
        tool="baleen",
    log:
        "logs/post_baleen_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_baleen_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule post_differr_sampled:
    input:
        'results/differr/{native}_{control}-{sample_size}-{n}/differr.bed'
    output:
        "results/differr/{native}_{control}-{sample_size}-{n}/differr_results.tsv"
    params:
        tool="differr",
    log:
        "logs/post_differr_sampled_format/{native}_{control}_{sample_size}_{n}.log"
    benchmark:
        "benchmarks/{native}_{control}_{sample_size}_{n}.post_differr_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule post_differr:
    input:
        'results/differr/{native}_{control}/differr.bed'
    output:
        "results/differr/{native}_{control}/differr_results.tsv"
    params:
        tool="differr",
    log:
        "logs/post_differr_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_differr_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


