rule samtools_filter_mapped:
    input:
        bam="results/alignments/{sample}.bam",
        csi="results/alignments/{sample}.bam.csi"
    output:
        "results/alignments/{sample}_filtered.bam",
        "results/alignments/{sample}_filtered.bam.csi"
    log:
        "logs/readfiltering/{sample}.log"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config['params']['samtools_filtering']
    shell:
        "samtools view -Sbh {params.extra} --write-index -o {output} {input.bam} 2>{log}"