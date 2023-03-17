rule slow5tools_fast5toslow5:
    input:
        fast5="data/{sample}/fast5/workspace"
    output:
        directory("results/slow5/{sample}")
    log:
        "logs/fast5toslow5/{sample}.log"
    threads: config['threads']['slow5tools']
    # conda:
    #     "../envs/slow5tools.yaml"
    container:
        "docker://btrspg/slow5tools:dev"
    shell:
        "slow5tools f2s {input.fast5} -d {output} -p {threads} "
        "2>{log}"

rule slow5tools_slow5toblow5:
    input:
        slow5="results/slow5/{sample}"
    output:
        "results/blow5/{sample}.blow5"
    log:
        "logs/slow5toblow5/{sample}.log"
    threads: config['threads']['slow5tools']
    # conda:
    #     "../envs/slow5tools.yaml"
    container:
        "docker://btrspg/slow5tools:dev"
    shell:
        "slow5tools merge {input.slow5} -o {output} -t{threads} "
        "2>{log} && rm -rf {input.slow5}"