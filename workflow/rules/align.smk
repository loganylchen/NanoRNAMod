rule minimap2_align:
    input:
        reference=config['reference']['transcriptome_fasta'],
        fastq="data/{sample}/fastq/pass.fq.gz"
    output:
        bam="results/alignments/{sample}.bam",
        csi="results/alignments/{sample}.bam.csi"
    log:
        "logs/minimap2_map/{sample}.log"
    benchmark:
        "benchmarks/{sample}.minimap2.benchmark.txt"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config["params"]["minimap2"]
    threads: config['threads']['minimap2']
    shell:
        "minimap2 -t {threads} {params.extra} {input.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} "
