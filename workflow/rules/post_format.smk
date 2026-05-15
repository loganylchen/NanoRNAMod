rule post_tandemmod:
    input:
        predictions="{project}/results/dataprep/{sample}_tandemmod_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_tandemmod_dataprep"
    output:
        "{project}/results/modifications/tandemmod/{sample}/tandemmod_results.tsv",
    params:
        tool="tandemmod",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_tandemmod_sampled_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_tandemmod_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule post_directrm:
    input:
        predictions="{project}/results/dataprep/{sample}_directrm_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_directrm_dataprep"
    output:
        "{project}/results/modifications/directrm/{sample}/directrm_results.tsv",
    params:
        tool="directrm",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_directrm_sampled_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_directrm_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule post_m6atm:
    input:
        predictions="{project}/results/dataprep/{sample}_m6atm_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_m6atm_dataprep"
    output:
        "{project}/results/modifications/m6atm/{sample}/m6atm_results.tsv",
    params:
        tool="m6atm",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_m6atm_sampled_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_m6atm_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule post_rnano:
    input:
        predictions="{project}/results/dataprep/{sample}_rnano_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_rnano_dataprep"
    output:
        "{project}/results/modifications/rnano/{sample}/rnano_results.tsv",
    params:
        tool="rnano",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_rnano_sampled_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_rnano_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


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


rule post_baleen:
    input:
        "{project}/results/baleen/{native}_{control}/transcript_mod_results.csv",
    output:
        "{project}/results/modifications/baleen/{native}_{control}/baleen_results.tsv",
    params:
        tool="baleen",
    container:
        get_container("python3")
    priority: 20
    log:
        "logs/{project}/post_baleen_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_baleen_sampled_format.benchmark.txt"
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


rule post_psipore:
    input:
        "{project}/results/psipore/{native}_{control}/psipore_results.tsv",
    output:
        "{project}/results/modifications/psipore/{native}_{control}/psipore_results.tsv",
    params:
        tool="psipore",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_psipore_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_psipore_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule post_nanopsu:
    input:
        predictions="{project}/results/dataprep/{sample}_nanopsu_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_nanopsu_dataprep",
    output:
        "{project}/results/modifications/nanopsu/{sample}/nanopsu_results.tsv",
    params:
        tool="nanopsu",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_nanopsu_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_nanopsu_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule post_nanomud:
    input:
        predictions="{project}/results/dataprep/{sample}_nanomud_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_nanomud_dataprep",
    output:
        "{project}/results/modifications/nanomud/{sample}/nanomud_results.tsv",
    params:
        tool="nanomud",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_nanomud_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_nanomud_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule post_penguin:
    input:
        predictions="{project}/results/dataprep/{sample}_penguin_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_penguin_dataprep",
    output:
        "{project}/results/modifications/penguin/{sample}/penguin_results.tsv",
    params:
        tool="penguin",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_penguin_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_penguin_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("default", 1)
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
