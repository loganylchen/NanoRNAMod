rule f5c_index:
    input:
        blow5="results/blow5/{sample}.blow5",
        fastq="results/fastq/{sample}.fq.gz",
    output:
        blow5_index="results/blow5/{sample}.blow5.idx",
        fastq_index=multiext("results/fastq/{sample}.fq.gz",'.index','.index.fai','.index.gzi','.index.readdb')
    log:
        "logs/f5c_index/{sample}.log"
    benchmark:
        "benchmarks/{sample}.f5c_index.benchmark.txt"
    conda:
        "../envs/f5c.yaml"
    threads: config['threads']['f5c']
    shell:
        "f5c index "
        "--slow5 {input.blow5} "
        "{input.fastq} "
        "-t {threads} "
        "2>{log} && "
        'echo -e "*\t{input.blow5}" > {input.fastq}.index.readdb'
