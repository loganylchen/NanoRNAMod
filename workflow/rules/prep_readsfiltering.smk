rule samtools_filter_mapped:
    input:
        bam="results/alignments/{sample}.bam",
        csi="results/alignments/{sample}.bam.csi"
    output:
        bam="results/alignments/{sample}_filtered.bam",
        bai="results/alignments/{sample}_filtered.bam.bai",
        csi="results/alignments/{sample}_filtered.bam.csi",
    log:
        "logs/readfiltering/{sample}.log"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config['params']['samtools_filtering']
    shell:
        "samtools view -Sbh {params.extra} --write-index -o {output.bam} {input.bam} 2>{log} && "
        "samtools index {output.bam} "

rule samtools_depth:
    input:
        bams = expand("results/alignments/{sample}_filtered.bam",sample=list(samples.index)),
        bai= expand("results/alignments/{sample}_filtered.bam.bai",sample=list(samples.index)),
    output:
        temp("results/bams.depth")
    log:
        "logs/depth.log"
    conda:
        "../envs/minimap2.yaml"
    shell:
        "samtools depth  -J {input.bams} > {output} 2>{log}  "

rule make_depth_table:
    input:
        "results/bams.depth"
    output:
        "results/depth_table.tsv"
    params:
        samples=list(samples.index)
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/depth.py"

rule samtools_filter_mapped_epi:
    input:
        bam="results/alignments/{sample}_3.2.4.bam",
        csi="results/alignments/{sample}_3.2.4.bam.csi"
    output:
        bam="results/alignments/{sample}_of.bam",
        bai="results/alignments/{sample}_of.bam.bai",
        csi="results/alignments/{sample}_of.bam.csi",
    log:
        "logs/readfiltering_3.2.4/{sample}.log"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config['params']['samtools_filtering']
    shell:
        "samtools view -Sbh {params.extra} --write-index -o {output.bam} {input.bam} 2>{log} && "
        "samtools index {output.bam} "


