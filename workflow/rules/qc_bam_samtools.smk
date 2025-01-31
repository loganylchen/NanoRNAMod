
rule qc_samtools:
    input:
        bam1="{project}/results/alignments/{sample}.splice.bam",
        bai1="{project}/results/alignments/{sample}.splice.bam.bai",
        bam2="{project}/results/alignments/{sample}.bam",
        bai2="{project}/results/alignments/{sample}.bam.bai",
    output:
        report1="{project}/results/qc/{sample}/{sample}_samtools_transcriptome_stats.txt",
        report2="{project}/results/qc/{sample}/{sample}_samtools_genome_stats.txt",
    benchmark:
        "benchmarks/{project}/{sample}.samtools_rnaseq.benchmark.txt"
    log:
        "logs/{project}/qc/{sample}_samtools.log",
    conda:
        "../envs/minimap2.yaml"
    priority: 50
    threads: config["threads"]["minimap2"]
    resources:
        mem_mb=1024 * 20,
    shell:
        "samtools stats -@ {threads} "
        "{input.bam1} > {output.report1} 2>>{log} && "
        "samtools stats -@ {threads} "
        "{input.bam2} > {output.report2} 2>>{log} "
