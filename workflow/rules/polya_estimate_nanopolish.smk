rule nanopolish_transcriptome_polya:
    input:
        bam="{project}/results/alignments/{sample}.bam",
        bai="{project}/results/alignments/{sample}.bam.bai",
        blow5_index="{project}/results/blow5/{sample}.blow5.idx",
        fastq_index=multiext(
            "{project}/results/fastq/{sample}.fq.gz",
            ".index",
            ".index.fai",
            ".index.gzi",
            ".index.readdb",
        ),
        fastq="{project}/results/fastq/{sample}.fq.gz",
    output:
        polya_estimate=KEEP_OR_NOT("{project}/results/polya/{sample}.transcriptome.nanopolish.polya.tsv.gz"),
    container:
        get_container("nanopolish")
    log:
        "logs/{project}/nanopolish_transcriptome_polya/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.nanopolish_transcriptome_polya.benchmark.txt"
    params:
        reference=config["reference"]["transcriptome_fasta"],
    resources:
        mem_mb = 1024 * 30
    threads: get_threads("nanopolish", 4)
    shell:
        "nanopolish polya "
        "--reads {input.fastq} "
        "--bam {input.bam} "
        "--genome {params.reference} -t {threads} 2>> {log} |gzip -c > {output.polya_estimate} "


rule nanopolish_genome_polya:
    input:
        bam="{project}/results/alignments/{sample}.splice.bam",
        bai="{project}/results/alignments/{sample}.splice.bam.bai",
        blow5_index="{project}/results/blow5/{sample}.blow5.idx",
        fastq_index=multiext(
            "{project}/results/fastq/{sample}.fq.gz",
            ".index",
            ".index.fai",
            ".index.gzi",
            ".index.readdb",
        ),
        fastq="{project}/results/fastq/{sample}.fq.gz",
    output:
        polya_estimate=KEEP_OR_NOT("{project}/results/polya/{sample}.genome.nanopolish.polya.tsv.gz"),
    container:
        get_container("nanopolish")
    log:
        "logs/{project}/nanopolish_genome_polya/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.nanopolish_genome_polya.benchmark.txt"
    params:
        reference=config["reference"]["genome_fasta"],
    resources:
        mem_mb = 1024 * 30
    threads: get_threads("nanopolish", 4)
    shell:
        "nanopolish polya "
        "--reads {input.fastq} "
        "--bam {input.bam} "
        "--genome {params.reference} -t {threads} 2>> {log} |gzip -c > {output.polya_estimate} "

