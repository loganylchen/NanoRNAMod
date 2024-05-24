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
        bams=["results/alignments/{sample}_filtered_"+f"{sample_size}_{n}.bam" for sample_size in config['sample_size'] for n in iter_number],
        bais=["results/alignments/{sample}_filtered_"+f"{sample_size}_{n}.bam.bai" for sample_size in config['sample_size'] for n in iter_number],
        read_names = ["results/alignments/{sample}_filtered_"+f"{sample_size}_{n}.txt" for sample_size in config['sample_size'] for n in iter_number],
    log:
        "logs/samplereads/{sample}.log"
    conda:
        "../envs/pysam.yaml"
    threads: config['threads']['sample_reads']
    params:
        sample_size=config['sample_size'],
        iter_number=iter_number
    script:
        "../scripts/sample_reads.py"
