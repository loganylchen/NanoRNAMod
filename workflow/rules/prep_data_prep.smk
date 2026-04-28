localrules:
    link_fastq,
    link_blow5,


rule link_fastq:
    input:
        fastq=get_raw_fastq,
    output:
        fastq="{project}/results/fastq/{sample}.fq.gz",
    params:
        relative_path=lambda w: "../../../" + get_raw_fastq(w),
    container:
        get_container("default")
    log:
        "logs/{project}/link_fastq/{sample}.log",
    threads: get_threads("default", 1)
    resources:
        mem_mb=1024 * 10,
    priority: 100
    shell:
        "ln -s {params.relative_path} {output.fastq} && "
        "echo `date` > {log} "


# rule link_fastq_old:
#     input:
#         old_fastq="data/{sample}/fastq/old.fq.gz",
#     output:
#         old_fastq="results/fastq/{sample}_3.2.4.fq.gz",
#     params:
#         old_relative_path="../../data/{sample}/fastq/old.fq.gz"
#     log:
#         "logs/link_fastq_old/{sample}.log"
#     shell:
#         "ln -s {params.old_relative_path} {output.old_fastq} && "
#         "echo `date` > {log} "


rule link_blow5:
    input:
        blow5=get_raw_blow5,
    output:
        blow5="{project}/results/blow5/{sample}.blow5",
    params:
        relative_path=lambda w: "../../../"+ get_raw_blow5(w),
    container:
        get_container("default")
    log:
        "logs/{project}/link_blow5/{sample}.log",
    threads: get_threads("default", 1)
    resources:
        mem_mb=1024 * 10,
    priority: 100
    shell:
        "ln -s {params.relative_path} {output.blow5} && "
        "echo `date` > {log} "
