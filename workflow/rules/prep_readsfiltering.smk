rule samtools_filter_mapped_epi:
    input:
        bam="{project}/results/alignments/{sample}_3.2.4.bam",
        csi="{project}/results/alignments/{sample}_3.2.4.bam.csi",
    output:
        bam=temp("{project}/results/alignments/{sample}_of.bam"),
        bai=temp("{project}/results/alignments/{sample}_of.bam.bai"),
        csi=temp("{project}/results/alignments/{sample}_of.bam.csi"),
    container:
        get_container("samtools")
    log:
        "logs/{project}/readfiltering_3.2.4/{sample}.log",
    params:
        extra=config["params"]["samtools_filtering"],
    shell:
        "samtools view -Sbh {params.extra} --write-index -o {output.bam} {input.bam} 2>{log} && "
        "samtools index {output.bam} "
