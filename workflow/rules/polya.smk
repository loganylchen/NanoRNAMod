rule nanopolish_polya:
    input:
        fastq="results/fastq/{sample}.fq.gz",
        bam="results/alignments/{sample}.splice.bam",
        index="results/fastq/{sample}.fq.gz.index"
    output:
        polya_tsv="results/polya/{sample}.tsv.gz",
    params:
        genome=config['reference']['genome_fasta']
    log:
        "logs/nanopolish_polya/{sample}.log"
    threads: config['threads']['nanopolish']
    conda:
        "../envs/nanopolish.yaml"
    shell:
        "nanopolish polya -t {threads} --reads {input.fastq} --bam {input.bam} --genome {params.genome} | "
        "gzip -c > {output.polya_tsv}  2>{log}"