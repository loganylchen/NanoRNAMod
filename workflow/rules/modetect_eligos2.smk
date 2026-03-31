rule eligos2:
    input:
        control_bam="{project}/results/alignments/{control}.bam",
        control_bai="{project}/results/alignments/{control}.bam.bai",
        native_bam="{project}/results/alignments/{native}.bam",
        native_bai="{project}/results/alignments/{native}.bam.bai",
        reference=config["reference"]["transcriptome_fasta"],
        region="{project}/results/eligos2/{native}_{control}.bed",
    output:
        directory=temp(directory("{project}/results/eligos2/{native}_{control}")),
        result=temp("{project}/results/eligos2/{native}_{control}/{native}_vs_{control}_on_{native}_{control}_baseExt0.txt"),
    params:
        prefix="{native}_{control}",
        extra=config["params"]["eligos2"],
    log:
        "logs/{project}/eligos2/{native}_{control}.log",
    threads: get_threads("eligos2", 4)
    benchmark:
        "benchmarks/{project}/{native}_{control}.eligos2.benchmark.txt"
    container:
        get_container("eligos2")
    resources:
        mem_mb = 1024 * 50
    priority: 10
    shell:
        "eligos2 pair_diff_mod "
        "-tbam {input.native_bam} "
        "-cbam {input.control_bam} "
        "-ref {input.reference} "
        "-t {threads} "
        "-reg {input.region} "
        "-o {output.directory} {params.extra} 1>{log} 2>&1"


rule eligos2_prep:
    input:
        control_bam="{project}/results/alignments/{control}.bam",
        control_bai="{project}/results/alignments/{control}.bam.bai",
        native_bam="{project}/results/alignments/{native}.bam",
        native_bai="{project}/results/alignments/{native}.bam.bai",
    output:
        region=temp("{project}/results/eligos2/{native}_{control}.bed"),
    container:
        get_container("bedtools")
    threads: get_threads("bedtools", 1)
    resources:
        mem_mb = 1024 * 50
    priority: 10
    shell:
        "bedtools bamtobed -bed12 -i {input.control_bam}  > {output.region}.tmp && "
        "bedtools bamtobed -bed12 -i {input.native_bam} >> {output.region}.tmp && "
        "bedtools sort -i {output.region}.tmp > {output.region}.tmp2 && "
        "bedtools merge -i {output.region}.tmp2 > {output.region} && "
        "rm {output.region}.tmp {output.region}.tmp2"
