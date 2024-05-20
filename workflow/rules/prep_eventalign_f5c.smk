rule f5c_eventalign_with_readname:
    input:
        fastq="results/fastq/{sample}.fq.gz",
        fastq_index=multiext("results/fastq/{sample}.fq.gz",'.index','.index.fai','.index.gzi','.index.readdb'),
        bam="results/alignments/{sample}_filtered.bam",
        bai="results/alignments/{sample}_filtered.bam.bai",
        index="results/blow5/{sample}.blow5.idx",
        blow5="results/blow5/{sample}.blow5",
    output:
        outfile=temp("results/eventalign/{sample}_full.tsv.bz2"),
        completion=temp("results/eventalign/{sample}_full.tsv.completed"),
    log:
        "logs/eventalign/{sample}_full.log"
    params:
        extra=config['params']['f5c_eventalign_full'],
        reference=config['reference']['transcriptome_fasta']
    benchmark:
        "benchmarks/{sample}.eventalign_full.benchmark.txt"
    conda:
        "../envs/f5c.yaml"
    threads: config['threads']['f5c']
    shell:
        "f5c eventalign -r {input.fastq} "
        "{params.extra} "
        "-b {input.bam} "
        "-g {params.reference} "
        "-t {threads} "
        "--slow5 {input.blow5} "
        "2>{log} | bzip2 -cz > {output.outfile}  && echo `date` > {output.completion} "

rule f5c_eventalign_full_sampled:
    input:
        fastq="results/fastq/{sample}.fq.gz",
        fastq_index=multiext("results/fastq/{sample}.fq.gz",'.index','.index.fai','.index.gzi','.index.readdb'),
        bam="results/alignments/{sample}_filtered_{sample_size}_{n}.bam",
        bai="results/alignments/{sample}_filtered_{sample_size}_{n}.bam.bai",
        index="results/blow5/{sample}.blow5.idx",
        blow5="results/blow5/{sample}.blow5",
    output:
        outfile=temp("results/eventalign/{sample}_full_{sample_size}_{n}.tsv.bz2"),
        completion=temp("results/eventalign/{sample}_full_{sample_size}_{n}.tsv.completed"),
    log:
        "logs/eventalign/{sample}_full_{sample_size}_{n}.log"
    params:
        extra=config['params']['f5c_eventalign_full'],
        reference=config['reference']['transcriptome_fasta']
    benchmark:
        "benchmarks/{sample}.eventalign_full_{sample_size}_{n}.benchmark.txt"
    conda:
        "../envs/f5c.yaml"
    threads: config['threads']['f5c']
    shell:
        "f5c eventalign -r {input.fastq} "
        "{params.extra} "
        "-b {input.bam} "
        "-g {params.reference} "
        "-t {threads} "
        "--slow5 {input.blow5} "
        "2>{log} | bzip2 -cz > {output.outfile}  && echo `date` > {output.completion} "


rule f5c_eventalign_simple:
    input:
        fastq="results/fastq/{sample}.fq.gz",
        fastq_index=multiext("results/fastq/{sample}.fq.gz",'.index','.index.fai','.index.gzi','.index.readdb'),
        bam="results/alignments/{sample}_filtered.bam",
        bai="results/alignments/{sample}_filtered.bam.bai",
        index="results/blow5/{sample}.blow5.idx",
        blow5="results/blow5/{sample}.blow5",

    output:
        outfile=temp("results/eventalign/{sample}_simple.tsv.bz2"),
        completion=temp("results/eventalign/{sample}_simple.tsv.completed"),
    log:
        "logs/eventalign/{sample}_xpore.log"
    params:
        extra=config['params']['f5c_eventalign_simple'],
        reference=config['reference']['transcriptome_fasta']
    benchmark:
        "benchmarks/{sample}.eventalign_simple.benchmark.txt"
    conda:
        "../envs/f5c.yaml"
    threads: config['threads']['f5c']
    shell:
        "f5c eventalign -r {input.fastq} "
        "-b {input.bam} "
        "-g {params.reference} "
        "{params.extra} "
        "-t {threads} "
        "--slow5 {input.blow5} "
        " 2>{log} | bzip2 -cz > {output.outfile} && echo `date` > {output.completion}"

rule f5c_eventalign_simple_sampled:
    input:
        fastq="results/fastq/{sample}.fq.gz",
        fastq_index=multiext("results/fastq/{sample}.fq.gz",'.index','.index.fai','.index.gzi','.index.readdb'),
        bam="results/alignments/{sample}_filtered_{sample_size}_{n}.bam",
        bai="results/alignments/{sample}_filtered_{sample_size}_{n}.bam.bai",
        index="results/blow5/{sample}.blow5.idx",
        blow5="results/blow5/{sample}.blow5",

    output:
        outfile=temp("results/eventalign/{sample}_simple_{sample_size}_{n}.tsv.bz2"),
        completion=temp("results/eventalign/{sample}_simple_{sample_size}_{n}.tsv.completed"),
    log:
        "logs/eventalign/{sample}_simple_{sample_size}_{n}.log"
    params:
        extra=config['params']['f5c_eventalign_simple'],
        reference=config['reference']['transcriptome_fasta']
    benchmark:
        "benchmarks/{sample}.eventalign_simple_{sample_size}_{n}.benchmark.txt"
    conda:
        "../envs/f5c.yaml"
    threads: config['threads']['f5c']
    shell:
        "f5c eventalign -r {input.fastq} "
        "-b {input.bam} "
        "-g {params.reference} "
        "{params.extra} "
        "-t {threads} "
        "--slow5 {input.blow5} "
        " 2>{log} | bzip2 -cz > {output.outfile} && echo `date` > {output.completion}"


