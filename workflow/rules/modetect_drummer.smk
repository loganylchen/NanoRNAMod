rule prep_drummer_region:
    input:
        control_bam = "results/alignments/{control}_filtered.bam",
        control_bai = "results/alignments/{control}_filtered.bam.bai",
        native_bam = "results/alignments/{native}_filtered.bam",
        native_bai = "results/alignments/{native}_filtered.bam.bai",
    output:
        region=temp("results/drummer/{native}_{control}/regions.bed")
    log:
        "logs/drummer/{native}_{control}/regions.log",
    threads: 1
    container:
        "docker://btrspg/drummer:latest"
    shell:
        "bedtools bamtobed -bed12 -i {input.control_bam} | cut -f1 | sort > {output.region}.tmp && "
        "bedtools bamtobed -bed12 -i {input.native_bam} | cut -f1 | sort >> {output.region}.tmp && "  
        "sort {output.region}.tmp | uniq -d > {output.region} && "
        "rm {output.region}.tmp"

rule prep_drummer_region_sampled:
    input:
        control_bam="results/alignments/{control}_filtered_{sample_size}_{n}.bam",
        control_bai="results/alignments/{control}_filtered_{sample_size}_{n}.bam.bai",
        native_bam="results/alignments/{native}_filtered_{sample_size}_{n}.bam",
        native_bai="results/alignments/{native}_filtered_{sample_size}_{n}.bam.bai",
    output:
        region=temp("results/drummer/{native}_{control}-{sample_size}-{n}/regions.bed")
    log:
        "logs/drummer/{native}_{control}-{sample_size}-{n}/regions.log"
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
        region="results/drummer/{native}_{control}/regions.bed",
        reference=config['reference']['transcriptome_fasta']
    output:
        directory("results/drummer/{native}_{control}/")
    log:
        stdout="logs/differr/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.differr.benchmark.txt"
    container:
        "docker://btrspg/drummer:latest"
    shell:
        "DEUMMER "
        "-r {input.reference} "
        "-l {input.region} "
        " -c {input.control_bam} "
        " -t {input.native_bam} "
        " -o {output} "
        " -a isoform  2>{log.stdout}"




rule drummer_sampled:
    input:
        control_bam="results/alignments/{control}_filtered_{sample_size}_{n}.bam",
        control_bai="results/alignments/{control}_filtered_{sample_size}_{n}.bam.bai",
        native_bam="results/alignments/{native}_filtered_{sample_size}_{n}.bam",
        native_bai="results/alignments/{native}_filtered_{sample_size}_{n}.bam.bai",
        reference=config['reference']['transcriptome_fasta']
    output:
        directory("results/drummer/{native}_{control}-{sample_size}-{n}/")
    log:
        stdout="logs/differr/{native}_{control}-{sample_size}-{n}.log"
    benchmark:
        "benchmarks/{native}_{control}-{sample_size}-{n}.differr.benchmark.txt"
    container:
        "docker://btrspg/drummer:latest"
    shell:
        "DEUMMER "
        "-r {input.reference} "
        "-l {input.region} "
        " -c {input.control_bam} "
        " -t {input.native_bam} "
        " -o {output} "
        " -a isoform  2>{log.stdout}"
