rule post_xpore:
    input:
        "{project}/results/xpore/{native}_{control}/diffmod.table",
    output:
        "{project}/results/modifications/xpore/{native}_{control}/xpore_results.tsv",
    params:
        tool="xpore",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_xpore_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_xpore_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 95
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule post_nanocompore:
    input:
        "{project}/results/nanocompore/{native}_{control}/nanocompore_results.tsv",
        "{project}/results/nanocompore/{native}_{control}",
    output:
        "{project}/results/modifications/nanocompore/{native}_{control}/nanocompore_results.tsv",
    params:
        tool="nanocompore",
    container:
        get_container("python3")
    priority: 86
    log:
        "logs/{project}/post_nanocompore_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_nanocompore_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule post_differr:
    input:
        "{project}/results/differr/{native}_{control}/differr.bed",
    output:
        "{project}/results/modifications/differr/{native}_{control}/differr_results.tsv",
    params:
        tool="differr",
    container:
        get_container("python3")
    priority: 75
    log:
        "logs/{project}/post_differr_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_differr_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule post_epinano:
    input:
        "{project}/results/epinano/{native}_{control}/epinano.delta-sum_err.prediction.csv",
        "{project}/results/epinano/{native}_{control}",
    output:
        "{project}/results/modifications/epinano/{native}_{control}/epinano_results.tsv",
    params:
        tool="epinano",
    container:
        get_container("python3")
    resources:
        mem_mb=1024 * 50,
    priority: 68
    threads: get_threads("python3", 1)
    log:
        "logs/{project}/post_epinano_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_epinano_sampled_format.benchmark.txt"
    script:
        "../scripts/format.py"


rule post_eligos2:
    input:
        "{project}/results/eligos2/{native}_{control}/{native}_vs_{control}_on_{native}_{control}_baseExt0.txt",
        "{project}/results/eligos2/{native}_{control}",
    output:
        "{project}/results/modifications/eligos2/{native}_{control}/eligos2_results.tsv",
    params:
        tool="eligos2",
    container:
        get_container("python3")
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("python3", 1)
    priority: 72
    log:
        "logs/{project}/post_eligos2_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_eligos2_sampled_format.benchmark.txt"
    script:
        "../scripts/format.py"


rule post_drummer:
    input:
        "{project}/results/drummer/{native}_{control}/",
    output:
        "{project}/results/modifications/drummer/{native}_{control}/drummer_results.tsv",
    params:
        tool="drummer",
    container:
        get_container("python3")
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("python3", 1)
    priority: 78
    log:
        "logs/{project}/post_drummer_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_drummer_sampled_format.benchmark.txt"
    script:
        "../scripts/format.py"


# pybaleen post-processing
rule post_pybaleen:
    input:
        result="{project}/results/pybaleen/{native}_{control}/site_results.tsv",
    output:
        result="{project}/results/modifications/pybaleen/{native}_{control}/pybaleen_results.tsv",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_pybaleen_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_pybaleen_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 100
    threads: get_threads("python3", 1)
    script:
        "../scripts/pybaleen_postprocess.py"
