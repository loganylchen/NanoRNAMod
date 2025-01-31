
rule qc_samtools:
    input:
        bam="{project}/results/alignments/{sample}.splice.bam",
        bai="{project}/results/alignments/{sample}.splice.bam.bai",
        bam2="{project}/results/alignments/{sample}.bam",
        bai2="{project}/results/alignments/{sample}.bam.bai",
    output:
        report1="{project}/results/qc/{sample}/{sample}_samtools_transcriptome_stats.txt",
        report2="{project}/results/qc/{sample}/{sample}_samtools_genome_stats.txt",
    params:
        output_dir=directory("{project}/results/qc/{sample}"),
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
        "{intput.bam1} > {otuput.report1} 2>>{log} "
        "samtools stats -@ {threads} "
        "{intput.bam2} > {otuput.report2} 2>>{log} "
