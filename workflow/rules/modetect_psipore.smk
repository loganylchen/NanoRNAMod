rule psipore_prep:
    input:
        control_bam="{project}/results/alignments/{control}_filtered.bam",
        control_bai="{project}/results/alignments/{control}_filtered.bam.bai",
        native_bam="{project}/results/alignments/{native}_filtered.bam",
        native_bai="{project}/results/alignments/{native}_filtered.bam.bai",
    output:
        region=temp("{project}/results/psipore/{native}_{control}_regions.bed"),
    container:
        get_container("bedtools")
    threads: get_threads("psipore", 1)
    resources:
        mem_mb = 1024
    priority: 10
    shell:
        "bedtools bamtobed -bed12 -i {input.control_bam} | cut -f1 | sort > {output.region}.tmp && "
        "bedtools bamtobed -bed12 -i {input.native_bam} | cut -f1 | sort >> {output.region}.tmp && "
        "sort {output.region}.tmp | uniq -d > {output.region} && "
        "rm {output.region}.tmp"


rule psipore_run:
    input:
        control_bam="{project}/results/alignments/{control}_filtered.bam",
        control_bai="{project}/results/alignments/{control}_filtered.bam.bai",
        native_bam="{project}/results/alignments/{native}_filtered.bam",
        native_bai="{project}/results/alignments/{native}_filtered.bam.bai",
        reference=config["reference"]["transcriptome_fasta"],
        region="{project}/results/psipore/{native}_{control}_regions.bed",
    output:
        result="{project}/results/psipore/{native}_{control}/psipore_results.tsv",
    container:
        get_container("psipore")
    threads: get_threads("psipore", 4)
    resources:
        mem_mb = 1024 * 50
    priority: 10
    params:
        extra=config["params"].get("psipore", ""),
        outdir="{project}/results/psipore/{native}_{control}",
    log:
        "logs/{project}/psipore/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.psipore.benchmark.txt"
    conda:
        "../envs/psipore.yaml"
    shell:
        "psipore compare "
        "--native {input.native_bam} "
        "--control {input.control_bam} "
        "--reference {input.reference} "
        "--regions {input.region} "
        "--threads {threads} "
        "--output {output.result} "
        "{params.extra} "
        "2>{log}"
