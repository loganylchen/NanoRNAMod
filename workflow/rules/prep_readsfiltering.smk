rule samtools_filter_mapped:
    input:
        bam="{project}/results/alignments/{sample}.bam",
        csi="{project}/results/alignments/{sample}.bam.csi"
    output:
        bam="{project}/results/alignments/{sample}_filtered.bam",
        bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        csi="{project}/results/alignments/{sample}_filtered.bam.csi",
    log:
        "logs/{project}/readfiltering/{sample}.log"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config['params']['samtools_filtering']
    shell:
        "samtools view -Sbh {params.extra} --write-index -o {output.bam} {input.bam} 2>{log} && "
        "samtools index {output.bam} "

rule samtools_depth:
    input:
        bams = expand("{project}/results/alignments/{sample}_filtered.bam",sample=list(samples.index)),
        bai= expand("{project}/results/alignments/{sample}_filtered.bam.bai",sample=list(samples.index)),
    output:
        temp("{project}/results/bams.depth")
    log:
        "logs/{project}/depth.log"
    conda:
        "../envs/minimap2.yaml"
    shell:
        "samtools depth  -J {input.bams} > {output} 2>{log}  "

rule make_depth_table:
    input:
        "{project}/results/bams.depth"
    output:
        "{project}/results/depth_table.tsv"
    params:
        samples=list(samples.index)
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/depth.py"

rule samtools_filter_mapped_epi:
    input:
        bam="{project}/results/alignments/{sample}_3.2.4.bam",
        csi="{project}/results/alignments/{sample}_3.2.4.bam.csi"
    output:
        bam="{project}/results/alignments/{sample}_of.bam",
        bai="{project}/results/alignments/{sample}_of.bam.bai",
        csi="{project}/results/alignments/{sample}_of.bam.csi",
    log:
        "logs/{project}/readfiltering_3.2.4/{sample}.log"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config['params']['samtools_filtering']
    shell:
        "samtools view -Sbh {params.extra} --write-index -o {output.bam} {input.bam} 2>{log} && "
        "samtools index {output.bam} "


