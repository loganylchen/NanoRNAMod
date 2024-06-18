rule prep_drummer_region:
    input:
        control_bam = "results/alignments/{control}_filtered.bam",
        control_bai = "results/alignments/{control}_filtered.bam.bai",
        native_bam = "results/alignments/{native}_filtered.bam",
        native_bai = "results/alignments/{native}_filtered.bam.bai",
    output:
        region=temp("results/drummer/{native}_{control}_regions.bed")
    threads: 1
    container:
        "docker://btrspg/drummer:latest"
    shell:
        "bedtools bamtobed -bed12 -i {input.control_bam} | cut -f1 | sort > {output.region}.tmp && "
        "bedtools bamtobed -bed12 -i {input.native_bam} | cut -f1 | sort >> {output.region}.tmp && "  
        "sort {output.region}.tmp | uniq -d > {output.region} && "
        "rm {output.region}.tmp"


rule drummer:
    input:
        control_bam="results/alignments/{control}_filtered.bam",
        control_bai="results/alignments/{control}_filtered.bam.bai",
        native_bam="results/alignments/{native}_filtered.bam",
        native_bai="results/alignments/{native}_filtered.bam.bai",
        region="results/drummer/{native}_{control}_regions.bed",
    params:
        reference = config['reference']['transcriptome_fasta'],
        extra=config['params']['drummer']
    output:
        outdir=directory("results/drummer/{native}_{control}/"),
        # summary="results/drummer/{native}_{control}/{control}_filtered-{native}_filtered/summary.txt"
    benchmark:
        "benchmarks/{native}_{control}.drummer.benchmark.txt"
    container:
        "docker://btrspg/drummer:latest"
    script:
        "../scripts/drummer.py"





