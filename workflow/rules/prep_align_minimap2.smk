rule minimap2_transcriptome_align:
    input:
        fastq="{project}/results/fastq/{sample}.fq.gz",
    output:
        bam=temp("{project}/results/alignments/{sample}.bam"),
        csi=temp("{project}/results/alignments/{sample}.bam.csi"),
        bai=temp("{project}/results/alignments/{sample}.bam.bai"),
    log:
        "logs/{project}/minimap2_transcriptome_alignment/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.minimap2_transcriptome_alignment.benchmark.txt"
    container:
        get_container("minimap2")
    params:
        extra=config["params"]["minimap2_transcriptome"],
        reference=config["reference"]["transcriptome_fasta"],
    resources:
        mem_mb = 1024 * 30
    threads: get_threads("minimap2", 4)
    shell:
        "minimap2 -t {threads} {params.extra} {params.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} && "
        "samtools index {output.bam}"


rule minimap2_genome_align:
    input:
        fastq="{project}/results/fastq/{sample}.fq.gz",
    output:
        bam=temp("{project}/results/alignments/{sample}.splice.bam"),
        csi=temp("{project}/results/alignments/{sample}.splice.bam.csi"),
        bai=temp("{project}/results/alignments/{sample}.splice.bam.bai"),
    log:
        "logs/{project}/minimap2_genome_alignment/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.minimap2_genome_alignment.benchmark.txt"
    container:
        get_container("minimap2")
    params:
        extra=config["params"]["minimap2_genome"],
        reference=config["reference"]["genome_fasta"],
    resources:
        mem_mb = 1024 * 30
    threads: get_threads("minimap2", 4)
    shell:
        "minimap2 -t {threads} {params.extra} {params.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} && "
        "samtools index {output.bam}"




rule minimap2_transcriptome_align_epi:
    input:
        fastq="{project}/results/fastq/{sample}_3.2.4.fq.gz",
    output:
        bam=temp("{project}/results/alignments/{sample}_3.2.4.bam"),
        csi=temp("{project}/results/alignments/{sample}_3.2.4.bam.csi"),
        bai=temp("{project}/results/alignments/{sample}_3.2.4.bam.bai"),
    log:
        "logs/{project}/minimap2_transcriptome_alignment_3.2.4/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.minimap2_transcriptome_alignment_3.2.4.benchmark.txt"
    container:
        get_container("minimap2")
    params:
        extra=config["params"]["minimap2_transcriptome"],
        reference=config["reference"]["transcriptome_fasta"],
    threads: get_threads("minimap2", 4)
    resources:
        mem_mb = 1024 * 30
    shell:
        "minimap2 -t {threads} {params.extra} {params.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} && "
        "samtools index {output.bam}"
