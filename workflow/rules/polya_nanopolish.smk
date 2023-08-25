rule nanopolish_polya:
    input:
        blow5_index="results/blow5/{sample}.blow5.idx",
        fastq="results/fastq/{sample}.fq.gz",
        fastq_index=multiext("results/fastq/{sample}.fq.gz",'.index','.index.fai','.index.gzi','.index.readdb'),
        bam="results/alignments/{sample}.splice.bam",
        bai="results/alignments/{sample}.splice.bam.bai",
    output:
        polya_tsv="results/polya/{sample}.tsv.gz",
    params:
        genome = config['reference']['genome_fasta']
    log:
        "logs/nanopolish_polya/{sample}.log"
    threads: config['threads']['nanopolish']
    conda:
        "../envs/nanopolish.yaml"
    shell:
        "nanopolish polya -t {threads} --reads {input.fastq} --bam {input.bam} --genome {params.genome} | "
        "gzip -c > {output.polya_tsv}  2>{log}"