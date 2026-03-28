rule prep_drummer_region:
    input:
        control_bam="{project}/results/alignments/{control}_filtered.bam",
        control_bai="{project}/results/alignments/{control}_filtered.bam.bai",
        native_bam="{project}/results/alignments/{native}_filtered.bam",
        native_bai="{project}/results/alignments/{native}_filtered.bam.bai",
    output:
        region=temp("{project}/results/drummer/{native}_{control}_regions.bed"),
    container:
        get_container("bedtools")
    threads: get_threads("drummer", 1)
    resources:
        mem_mb = 1024 * 10
    priority: 10
    shell:
        "bedtools bamtobed -bed12 -i {input.control_bam} | cut -f1 | sort > {output.region}.tmp && "
        "bedtools bamtobed -bed12 -i {input.native_bam} | cut -f1 | sort >> {output.region}.tmp && "
        "sort {output.region}.tmp | uniq -d > {output.region} && "
        "rm {output.region}.tmp"


rule drummer:
    input:
        control_bam="{project}/results/alignments/{control}_filtered.bam",
        control_bai="{project}/results/alignments/{control}_filtered.bam.bai",
        native_bam="{project}/results/alignments/{native}_filtered.bam",
        native_bai="{project}/results/alignments/{native}_filtered.bam.bai",
        region="{project}/results/drummer/{native}_{control}_regions.bed",
    params:
        reference=config["reference"]["transcriptome_fasta"],
        extra=config["params"]["drummer"],
    container:
        get_container("drummer")
    priority: 10
    log:
        "logs/{project}/drummer/{native}_{control}.log",
    resources:
        mem_mb = 1024 * 20
    output:
        outdir=KEEP_OR_NOT(directory("{project}/results/drummer/{native}_{control}/")),
        # summary="results/drummer/{native}_{control}/{control}_filtered-{native}_filtered/summary.txt"
    benchmark:
        "benchmarks/{project}/{native}_{control}.drummer.benchmark.txt"
    script:
        "../scripts/drummer.py"
