rule link_fastq:
    input:
        fastq="data/{sample}/fastq/pass.fq.gz"
    output:
        fastq="results/fastq/{sample}.fq.gz",
    params:
        relative_path="../../data/{sample}/fastq/pass.fq.gz"
    log:
        "logs/link_fastq/{sample}.log"
    shell:
        "ln -s {params.relative_path} {output.fastq} && "
        "echo `date` > {log} "