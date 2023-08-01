rule minimap2_transcriptome_align:
    input:
        reference=config['reference']['transcriptome_fasta'],
        fastq="data/{sample}/fastq/pass.fq.gz"
    output:
        bam="results/alignments/{sample}.bam",
        csi="results/alignments/{sample}.bam.csi"
    log:
        "logs/minimap2_map/{sample}.log"
    benchmark:
        "benchmarks/{sample}.minimap2.benchmark.txt"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config["params"]["minimap2_transcriptome"]
    threads: config['threads']['minimap2']
    shell:
        "minimap2 -t {threads} {params.extra} {input.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} "



rule minimap2_genome_align:
    input:
        reference=config['reference']['genome_fasta'],
        fastq="results/fastq/{sample}.fq.gz",
    output:
        bam="results/alignments/{sample}.splice.bam",
        csi="results/alignments/{sample}.splice.bam.csi"
    log:
        "logs/minimap2_map/{sample}.splice.log"
    benchmark:
        "benchmarks/{sample}.minimap2.splice.benchmark.txt"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config["params"]["minimap2_genome"]
    threads: config['threads']['minimap2']
    shell:
        "minimap2 -t {threads} {params.extra} {input.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} "
