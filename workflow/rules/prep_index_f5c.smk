rule f5c_index:
    input:
        blow5="{project}/results/blow5/{sample}.blow5",
        fastq="{project}/results/fastq/{sample}.fq.gz",
    output:
        blow5_index=temp("{project}/results/blow5/{sample}.blow5.idx"),
        fastq_index=temp(
            multiext(
                "{project}/results/fastq/{sample}.fq.gz",
                ".index",
                ".index.fai",
                ".index.gzi",
                ".index.readdb",
            )
        ),
    container:
        get_container("f5c")
    log:
        "logs/{project}/f5c_index/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.f5c_index.benchmark.txt"
    threads: get_threads("f5c", 4)
    resources:
        mem_mb = 1024 * 50,
    priority: 100
    shell:
        "f5c index "
        "--slow5 {input.blow5} "
        "{input.fastq} "
        "-t {threads} "
        "2>{log} && "
        'echo -e "*\t{input.blow5}" > {input.fastq}.index.readdb'
