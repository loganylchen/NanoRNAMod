rule f5c_index:
    input:
        blow5="results/blow5/{sample}.blow5",
        fastq="data/{sample}/fastq/pass.fq.gz",
    output:
        "results/blow5/{sample}.blow5.idx"
    log:
        "logs/f5c_index/{sample}.log"
    benchmark:
        "benchmarks/{sample}.f5cindex.benchmark.txt"
    conda:
        "../envs/f5c.yaml"
    threads: config['threads']['f5c']
    shell:
        "f5c index "
        "--slow5 {input.blow5} "
        "{input.fastq} "
        "-t {threads} "
        "2>{log}"


rule f5c_eventalign_nanocompore:
    input:
        fastq="data/{sample}/fastq/pass.fq.gz",
        bam="results/alignments/{sample}_filtered.bam",
        csi="results/alignments/{sample}_filtered.bam.csi",
        index="results/blow5/{sample}.blow5.idx",
        blow5="results/blow5/{sample}.blow5",
        reference=config['reference']['transcriptome_fasta']
    output:
        outfile="results/eventalign/{sample}_nanocompore.tsv.gz",
        completion="results/eventalign/{sample}_nanocompore.tsv.completed",
    log:
        "logs/eventalign/{sample}_nanocompore.log"
    params:
        extra=config['params']['f5c_eventalign_nanocompore']
    benchmark:
        "benchmarks/{sample}.eventalign_nanocompore.benchmark.txt"
    conda:
        "../envs/f5c.yaml"
    threads: config['threads']['f5c']
    shell:
        "f5c eventalign -r {input.fastq} "
        "{params.extra} "
        "-b {input.bam} "
        "-g {input.reference} "
        "-t {threads} "
        "--slow5 {input.blow5} "
        "2>{log} | gzip -c > {output.outfile}  && echo `date` > {output.completion} "


rule f5c_eventalign_xpore:
    input:
        fastq="data/{sample}/fastq/pass.fq.gz",
        bam="results/alignments/{sample}_filtered.bam",
        csi="results/alignments/{sample}_filtered.bam.csi",
        index="results/blow5/{sample}.blow5.idx",
        blow5="results/blow5/{sample}.blow5",
        reference=config['reference']['transcriptome_fasta']
    output:
        outfile="results/eventalign/{sample}_xpore.tsv.gz",
        completion="results/eventalign/{sample}_xpore.tsv.completed",
    log:
        "logs/eventalign/{sample}_xpore.log"
    params:
        extra=config['params']['f5c_eventalign_xpore']
    benchmark:
        "benchmarks/{sample}.eventalign_xpore.benchmark.txt"
    conda:
        "../envs/f5c.yaml"
    threads: config['threads']['f5c']
    shell:
        "f5c eventalign -r {input.fastq} "
        "-b {input.bam} "
        "-g {input.reference} "
        "{params.extra} "
        "-t {threads} "
        "--slow5 {input.blow5} "
        " 2>{log} | gzip -c > {output.outfile} && echo `date` > {output.completion}"