rule link_fastq:
    input:
        fastq="data/{sample}/fastq/pass.fq.gz",
        old_fastq="data/{sample}/fastq/old.fq.gz",
    output:
        fastq="results/fastq/{sample}.fq.gz",
        old_fastq="results/fastq/{sample}_3.2.4.fq.gz",
    params:
        relative_path="../../data/{sample}/fastq/pass.fq.gz",
        old_relative_path="../../data/{sample}/fastq/old.fq.gz"
    log:
        "logs/link_fastq/{sample}.log"
    shell:
        "ln -s {params.relative_path} {output.fastq} && "
        "ln -s {params.old_relative_path} {output.old_fastq} && "
        "echo `date` > {log} "



rule link_blow5:
    input:
        blow5="data/{sample}/blow5/nanopore.blow5"
    output:
        blow5="results/blow5/{sample}.blow5",
    params:
        relative_path="../../data/{sample}/blow5/nanopore.blow5"
    log:
        "logs/link_blow5/{sample}.log"
    shell:
        "ln -s {params.relative_path} {output.blow5} && "
        "echo `date` > {log} "







