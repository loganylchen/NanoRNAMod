rule f5c_index:
    input:
        blow5="{project}/results/blow5/{sample}.blow5",
        fastq="{project}/results/fastq/{sample}.fq.gz",
    output:
        blow5_index="{project}/results/blow5/{sample}.blow5.idx",
        fastq_index=multiext(
            "{project}/results/fastq/{sample}.fq.gz",
            ".index",
            ".index.fai",
            ".index.gzi",
            ".index.readdb",
        ),
    log:
        "logs/{project}/f5c_index/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.f5c_index.benchmark.txt"
    container:
        "docker://btrspg/f5c:dev"
    threads: config["threads"]["f5c"]
    resources:
        mem_mb = 1024 * 50,
    shell:
        "f5c index "
        "--slow5 {input.blow5} "
        "{input.fastq} "
        "-t {threads} "
        "2>{log} && "
        'echo -e "*\t{input.blow5}" > {input.fastq}.index.readdb'
