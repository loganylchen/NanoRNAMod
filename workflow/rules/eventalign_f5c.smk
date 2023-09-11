rule f5c_eventalign_nanocompore:
    input:
        fastq="data/{sample}/fastq/pass.fq.gz",
        bam="results/alignments/{sample}_filtered.bam",
        csi="results/alignments/{sample}_filtered.bam.csi",
        index="results/blow5/{sample}.blow5.idx",
        blow5="results/blow5/{sample}.blow5",

    output:
        outfile="results/eventalign/{sample}_nanocompore.tsv.gz",
        completion="results/eventalign/{sample}_nanocompore.tsv.completed",
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
        fastq="data/{sample}/fastq/pass.fq.gz",
        bam="results/alignments/{sample}_filtered.bam",
        csi="results/alignments/{sample}_filtered.bam.csi",
        index="results/blow5/{sample}.blow5.idx",
        blow5="results/blow5/{sample}.blow5",

    output:
        outfile="results/eventalign/{sample}_xpore.tsv.gz",
        completion="results/eventalign/{sample}_xpore.tsv.completed",
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
        fastq="data/{sample}/fastq/pass.fq.gz",
        bam="results/alignments/{sample}_filtered_sampling.bam",
        csi="results/alignments/{sample}_filtered_sampling.bam.bai",
        index="results/blow5/{sample}.blow5.idx",
        blow5="results/blow5/{sample}.blow5",
    output:
        outfile="results/eventalign/{sample}_baleen.tsv.bz2",
        completion="results/eventalign/{sample}_baleen.completed",
    log:
        "logs/eventalign/{sample}_baleen.log"
    params:
        extra=config['params']['f5c_eventalign_baleen'],
        reference=config['reference']['transcriptome_fasta']
    benchmark:
        "benchmarks/{sample}.eventalign_baleen.benchmark.txt"
    # container:
    #     "docker://btrspg/f5c:dev"
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
        # "--edparam {wildcards.w1},$(({wildcards.w1}+{wildcards.w2})),{wildcards.threshold},9.0,{wildcards.peak} "
        "2>{log} | bzip2 -cz > {output.outfile}  && echo `date` > {output.completion} "