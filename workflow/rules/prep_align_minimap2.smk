rule minimap2_transcriptome_align:
    input:
        fastq="results/fastq/{sample}.fq.gz",
    output:
        bam="results/alignments/{sample}.bam",
        csi="results/alignments/{sample}.bam.csi",
        bai="results/alignments/{sample}.bam.bai",
    log:
        "logs/minimap2_transcriptome_alignment/{sample}.log"
    benchmark:
        "benchmarks/{sample}.minimap2_transcriptome_alignment.benchmark.txt"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config["params"]["minimap2_transcriptome"],
        reference=config['reference']['transcriptome_fasta'],
    threads: config['threads']['minimap2']
    shell:
        "minimap2 -t {threads} {params.extra} {params.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} && "
        "samtools index {output.bam}"
