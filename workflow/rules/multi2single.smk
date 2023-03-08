rule multi_to_single:
    input:
        "data/{sample}/fast5/workspace"
    output:
        outdir=directory("results/multi_to_single/{sample}"),
        completion="results/multi_to_single/{sample}/m2s_complete.txt"
    log:
        "logs/multi_to_single/{sample}.log"
    benchmark:
        "benchmarks/{sample}.multi_to_single.benchmark.txt"
    threads: config['threads']['multi_to_single']
    conda:
        "../envs/ont-fast5-api.yaml"
    shell:
        "multi_to_single_fast5 "
        "-i {input} "
        "-s {output.outdir} "
        "-t {threads} --recursive 2>{log} && echo `date` > {output.completion}"