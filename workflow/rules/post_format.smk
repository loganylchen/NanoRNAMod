rule post_xpore_sampled:
    input:
        "results/xpore/{native}_{control}-{sample_size}-{n}/diffmod.table",
    output:
        "results/modifications/xpore/{native}_{control}-{sample_size}-{n}/xpore_results.tsv"
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
        "results/modifications/xpore/{native}_{control}/xpore_results.tsv"
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


rule post_nanocompore_sampled:
    input:
        "results/nanocompore/{native}_{control}-{sample_size}-{n}/nanocompore_results.tsv",
    output:
        "results/modifications/nanocompore/{native}_{control}-{sample_size}-{n}/nanocompore_results.tsv"
    params:
        tool="nanocompore",
    log:
        "logs/post_nanocompore_sampled_format/{native}_{control}_{sample_size}_{n}.log"
    benchmark:
        "benchmarks/{native}_{control}_{sample_size}_{n}.post_nanocompore_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule post_nanocompore:
    input:
        "results/nanocompore/{native}_{control}/nanocompore_results.tsv",
    output:
        "results/modifications/nanocompore/{native}_{control}/nanocompore_results.tsv"
    params:
        tool="nanocompore",
    log:
        "logs/post_nanocompore_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_nanocompore_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule post_baleen_sampled:
    input:
        'results/baleen/{native}_{control}-{sample_size}-{n}/transcripts.csv'
    output:
        "results/modifications/baleen/{native}_{control}-{sample_size}-{n}/baleen_results.tsv"
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
        "results/modifications/baleen/{native}_{control}/baleen_results.tsv"
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
        "results/modifications/differr/{native}_{control}-{sample_size}-{n}/differr_results.tsv"
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
        "results/modifications/differr/{native}_{control}/differr_results.tsv"
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


rule post_epinano:
    input:
        'results/epinano/{native}_{control}/epinano.delta-sum_err.prediction.csv'
    output:
        "results/modifications/epinano/{native}_{control}/epinano_results.tsv"
    params:
        tool="epinano",
    log:
        "logs/post_epinano_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_epinano_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule post_epinano_sampled:
    input:
        'results/epinano/{native}_{control}-{sample_size}-{n}/epinano.delta-sum_err.prediction.csv'
    output:
        "results/modifications/epinano/{native}_{control}-{sample_size}-{n}/epinano_results.tsv"
    params:
        tool="epinano",
    log:
        "logs/post_epinano_sampled_format/{native}_{control}-{sample_size}-{n}.log"
    benchmark:
        "benchmarks/{native}_{control}-{sample_size}-{n}.post_epinano_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule post_eligos2:
    input:
        "results/eligos2/{native}_{control}/{native}_filtered_vs_{control}_filtered_on_{native}_{control}_baseExt0.txt"
    output:
        "results/modifications/eligos2/{native}_{control}/eligos2_results.tsv"
    params:
        tool="eligos2",
    log:
        "logs/post_eligos2_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_eligos2_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule post_eligos2_sampled:
    input:
        "results/eligos2/{native}_{control}-{sample_size}-{n}/{native}_filtered_{sample_size}_{n}_vs_{control}_filtered_{sample_size}_{n}_on_{native}_filtered_{sample_size}_{n}_{control}_filtered_{sample_size}_{n}_baseExt0.txt",
    output:
        "results/modifications/eligos2/{native}_{control}-{sample_size}-{n}/eligos2_results.tsv"
    params:
        tool="eligos2",
    log:
        "logs/post_eligos2_sampled_format/{native}_{control}-{sample_size}-{n}.log"
    benchmark:
        "benchmarks/{native}_{control}-{sample_size}-{n}.post_eligos2_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule post_drummer:
    input:
        "results/drummer/{native}_{control}/"
    output:
        "results/modifications/drummer/{native}_{control}/drummer_results.tsv"
    params:
        tool="drummer",
    log:
        "logs/post_drummer_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_drummer_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule post_drummer_sampled:
    input:
        "results/drummer/{native}_{control}-{sample_size}-{n}/"
    output:
        "results/modifications/drummer/{native}_{control}-{sample_size}-{n}/drummer_results.tsv"
    params:
        tool="drummer",
    log:
        "logs/post_drummer_sampled_format/{native}_{control}-{sample_size}-{n}.log"
    benchmark:
        "benchmarks/{native}_{control}-{sample_size}-{n}.post_drummer_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"