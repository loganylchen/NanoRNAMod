rule f5c_eventalign_nanocompore:
    input:
        fastq="results/fastq/{sample}.fq.gz",
        fastq_index=multiext("results/fastq/{sample}.fq.gz",'.index','.index.fai','.index.gzi','.index.readdb'),
        bam="results/alignments/{sample}_filtered.bam",
        bai="results/alignments/{sample}_filtered.bam.bai",
        index="results/blow5/{sample}.blow5.idx",
        blow5="results/blow5/{sample}.blow5",

    output:
        outfile=temp("results/eventalign/{sample}_nanocompore.tsv.gz"),
        completion=temp("results/eventalign/{sample}_nanocompore.tsv.completed"),
    log:
        "logs/eventalign/{sample}_nanocompore.log"
    params:
        extra=config['params']['f5c_eventalign_nanocompore'],
        reference=config['reference']['transcriptome_fasta']
    benchmark:
        "benchmarks/{sample}.eventalign_nanocompore.benchmark.txt"
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
        "2>{log} | gzip -c > {output.outfile}  && echo `date` > {output.completion} "


rule f5c_eventalign_xpore:
    input:
        fastq="results/fastq/{sample}.fq.gz",
        fastq_index=multiext("results/fastq/{sample}.fq.gz",'.index','.index.fai','.index.gzi','.index.readdb'),
        bam="results/alignments/{sample}_filtered.bam",
        bai="results/alignments/{sample}_filtered.bam.bai",
        index="results/blow5/{sample}.blow5.idx",
        blow5="results/blow5/{sample}.blow5",

    output:
        outfile=temp("results/eventalign/{sample}_xpore.tsv.gz"),
        completion=temp("results/eventalign/{sample}_xpore.tsv.completed"),
    log:
        "logs/eventalign/{sample}_xpore.log"
    params:
        extra=config['params']['f5c_eventalign_xpore'],
        reference=config['reference']['transcriptome_fasta']
    benchmark:
        "benchmarks/{sample}.eventalign_xpore.benchmark.txt"
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
        " 2>{log} | gzip -c > {output.outfile} && echo `date` > {output.completion}"


rule f5c_eventalign_baleen:
    input:
        fastq="results/fastq/{sample}.fq.gz",
        fastq_index=multiext("results/fastq/{sample}.fq.gz",'.index','.index.fai','.index.gzi','.index.readdb'),
        bam="results/alignments/{sample}.realign.bam" if config['realign'] else "results/alignments/{sample}_filtered.bam",
        bai="results/alignments/{sample}.realign.bam.bai" if config['realign'] else "results/alignments/{sample}_filtered.bam.bai",
        index="results/blow5/{sample}.blow5.idx",
        blow5="results/blow5/{sample}.blow5",
    output:
        outfile="results/eventalign/{sample}_baleen.tsv.bz2",
        completion="results/eventalign/{sample}_baleen.tsv.completed",
    log:
        "logs/eventalign/{sample}_baleen.log"
    params:
        extra=config['params']['f5c_eventalign_baleen'],
        reference=config['reference']['transcriptome_fasta']
    benchmark:
        "benchmarks/{sample}.eventalign_baleen.benchmark.txt"
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
        "2>{log} | bzip2 -cz  > {output.outfile}  && echo `date` > {output.completion} "
