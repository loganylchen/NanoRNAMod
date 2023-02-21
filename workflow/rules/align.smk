rule align:
    input:
        reference=config['reference'],
        fastq="data/{sample}/fastq/pass.fq.gz"
    output:
        bam="results/mapped_reads/{sample}.bam"
    log:
        "logs/minimap2_map/{sample}.log"
    benchmark:
        "benchmarks/{sample}.minimap2.benchmark.txt"
    threads: 12
    shell:
        "minimap2 -t {threads} -ax map-ont -L --secondary=no -N 1 {input.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} "
