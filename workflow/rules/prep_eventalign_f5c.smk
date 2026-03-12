rule f5c_eventalign_with_readname:
    input:
        fastq="{project}/results/fastq/{sample}.fq.gz",
        fastq_index=multiext(
            "{project}/results/fastq/{sample}.fq.gz",
            ".index",
            ".index.fai",
            ".index.gzi",
            ".index.readdb",
        ),
        bam="{project}/results/alignments/{sample}_filtered.bam",
        bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        index="{project}/results/blow5/{sample}.blow5.idx",
        blow5="{project}/results/blow5/{sample}.blow5",
    output:
        outfile=KEEP_OR_NOT("{project}/results/eventalign/{sample}_full.tsv.bz2"),
        completion=KEEP_OR_NOT(
            "{project}/results/eventalign/{sample}_full.tsv.completed"
        ),
    container:
        get_container("f5c")
    priority: 20
    log:
        "logs/{project}/eventalign/{sample}_full.log",
    params:
        extra=config["params"]["f5c_eventalign_full"],
        reference=config["reference"]["transcriptome_fasta"],
    benchmark:
        "benchmarks/{project}/{sample}.eventalign_full.benchmark.txt"
    conda:
        "../envs/f5c.yaml"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("f5c", 4)
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
        fastq="{project}/results/fastq/{sample}.fq.gz",
        fastq_index=multiext(
            "{project}/results/fastq/{sample}.fq.gz",
            ".index",
            ".index.fai",
            ".index.gzi",
            ".index.readdb",
        ),
        bam="{project}/results/alignments/{sample}_filtered.bam",
        bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        index="{project}/results/blow5/{sample}.blow5.idx",
        blow5="{project}/results/blow5/{sample}.blow5",
    output:
        outfile=KEEP_OR_NOT("{project}/results/eventalign/{sample}_simple.tsv.bz2"),
        completion=KEEP_OR_NOT(
            "{project}/results/eventalign/{sample}_simple.tsv.completed"
        ),
    container:
        get_container("f5c")
    log:
        "logs/{project}/eventalign/{sample}_xpore.log",
    params:
        extra=config["params"]["f5c_eventalign_simple"],
        reference=config["reference"]["transcriptome_fasta"],
    benchmark:
        "benchmarks/{project}/{sample}.eventalign_simple.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    conda:
        "../envs/f5c.yaml"
    threads: get_threads("f5c", 4)
    shell:
        "f5c eventalign -r {input.fastq} "
        "-b {input.bam} "
        "-g {params.reference} "
        "{params.extra} "
        "-t {threads} "
        "--slow5 {input.blow5} "
        " 2>{log} | bzip2 -cz > {output.outfile} && echo `date` > {output.completion}"
