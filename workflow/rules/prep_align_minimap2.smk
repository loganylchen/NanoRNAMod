rule minimap2_transcriptome_align:
    input:
        fastq="{project}/results/fastq/{sample}.fq.gz",
    output:
        bam="{project}/results/alignments/{sample}.bam",
        csi="{project}/results/alignments/{sample}.bam.csi",
        bai="{project}/results/alignments/{sample}.bam.bai",
    log:
        "logs/{project}/minimap2_transcriptome_alignment/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.minimap2_transcriptome_alignment.benchmark.txt"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config["params"]["minimap2_transcriptome"],
        reference=config["reference"]["transcriptome_fasta"],
    threads: config["threads"]["minimap2"]
    shell:
        "minimap2 -t {threads} {params.extra} {params.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} && "
        "samtools index {output.bam}"


rule minimap2_transcriptome_align_epi:
    input:
        fastq="{project}/results/fastq/{sample}_3.2.4.fq.gz",
    output:
        bam="{project}/results/alignments/{sample}_3.2.4.bam",
        csi="{project}/results/alignments/{sample}_3.2.4.bam.csi",
        bai="{project}/results/alignments/{sample}_3.2.4.bam.bai",
    log:
        "logs/{project}/minimap2_transcriptome_alignment_3.2.4/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.minimap2_transcriptome_alignment_3.2.4.benchmark.txt"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config["params"]["minimap2_transcriptome"],
        reference=config["reference"]["transcriptome_fasta"],
    threads: config["threads"]["minimap2"]
    shell:
        "minimap2 -t {threads} {params.extra} {params.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} && "
        "samtools index {output.bam}"
