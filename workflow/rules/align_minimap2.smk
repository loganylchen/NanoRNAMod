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

rule minimap2_list_align:
    input:
        mapping_list="results/read_assignment/{sample}.list",
        mapping_dir="results/read_assignment/{sample}_tmp",
    output:
        bam_list="results/read_assignment/{sample}.bamlist",
        bam="results/alignments/{sample}.realign.bam",
        bai="results/alignments/{sample}.realign.bam.bai",
    log:
        "logs/minimap2_list_alignment/{sample}.log"
    benchmark:
        "benchmarks/{sample}.minimap2_list_alignment.benchmark.txt"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config["params"]["minimap2_transcriptome"],
        transcriptome_fasta=config['reference']['transcriptome_fasta']
    threads: config['threads']['minimap2']
    script:
        "../scripts/mapping_list.bash"


rule minimap2_genome_align:
    input:
        fastq="results/fastq/{sample}.fq.gz",
    output:
        bam="results/alignments/{sample}.splice.bam",
        csi="results/alignments/{sample}.splice.bam.csi",
        bai="results/alignments/{sample}.splice.bam.bai"
    log:
        "logs/minimap2_genome_alignment/{sample}.log"
    benchmark:
        "benchmarks/{sample}.minimap2_genome_alignment.benchmark.txt"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config["params"]["minimap2_genome"],
        reference=config['reference']['genome_fasta'],
    threads: config['threads']['minimap2']
    shell:
        "minimap2 -t {threads} {params.extra} {params.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} && "
        "samtools index {output.bam}"
