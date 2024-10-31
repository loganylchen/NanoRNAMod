rule differr:
    input:
        control_bam="{project}/results/alignments/{control}_filtered.bam",
        control_bai="{project}/results/alignments/{control}_filtered.bam.bai",
        native_bam="{project}/results/alignments/{native}_filtered.bam",
        native_bai="{project}/results/alignments/{native}_filtered.bam.bai",
        reference=config["reference"]["transcriptome_fasta"],
    output:
        KEEP_OR_NOT("{project}/results/differr/{native}_{control}/differr.bed"),
    params:
        prefix="",
        extra=config["params"]["differr"],
    log:
        stdout="logs/{project}/differr/{native}_{control}.log",
    threads: config["threads"]["differr"]
    benchmark:
        "benchmarks/{project}/{native}_{control}.differr.benchmark.txt"
    container:
        "docker://btrspg/differr:latest"
    shell:
        "differr "
        " -a {input.control_bam} "
        " -b {input.native_bam} "
        " -r {input.reference} "
        " -o {output} "
        " {params.extra} "
        " -p {threads} 2>{log}"
