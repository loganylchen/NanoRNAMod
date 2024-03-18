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

rule sample_reads:
    input:
        bam="results/alignments/{sample}_filtered.bam",
        bai="results/alignments/{sample}_filtered.bam.bai",
    output:
        bams=["results/alignments/{sample}_filtered_"+f"{sample_size}.bam" for sample_size in config['sample_size']],
        bais=["results/alignments/{sample}_filtered_"+f"{sample_size}.bam.bai" for sample_size in config['sample_size']],
        read_names = ["results/alignments/{sample}_filtered_"+f"{sample_size}.txt" for sample_size in config['sample_size']],
    log:
        "logs/samplereads/{sample}.log"
    conda:
        "../envs/pysam.yaml"
    threads: config['threads']['sample_reads']
    params:
        sample_size=config['sample_size']
    script:
        "../scripts/sample_reads.py"
