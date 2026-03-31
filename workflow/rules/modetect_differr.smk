rule differr:
    input:
        control_bam="{project}/results/alignments/{control}.bam",
        control_bai="{project}/results/alignments/{control}.bam.bai",
        native_bam="{project}/results/alignments/{native}.bam",
        native_bai="{project}/results/alignments/{native}.bam.bai",
        reference=config["reference"]["transcriptome_fasta"],
    output:
        temp("{project}/results/differr/{native}_{control}/differr.bed"),
    params:
        prefix="",
        extra=config["params"]["differr"],
    container:
        get_container("differr")
    log:
        "logs/{project}/differr/{native}_{control}.log",
    threads: get_threads("differr", 4)
    resources:
        mem_mb = 1024 * 50
    priority: 10
    benchmark:
        "benchmarks/{project}/{native}_{control}.differr.benchmark.txt"
    shell:
        "differr "
        " -a {input.control_bam} "
        " -b {input.native_bam} "
        " -r {input.reference} "
        " -o {output} "
        " {params.extra} "
        " -p {threads} 1>{log} 2>&1"
