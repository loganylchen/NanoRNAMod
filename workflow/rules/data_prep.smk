rule link_fastq:
    input:
        fastq="data/{sample}/fastq/pass.fq.gz"
    output:
        fastq="results/fastq/{sample}.fq.gz",
    log:
        "logs/cp/{sample}.log"
    shell:
        "ln -s {input.fastq} {output.fastq}"

