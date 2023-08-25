rule slow5tools_f2s:
    input:
        fast5="data/{sample}/fast5/workspace"
    output:
        outdir=temp(directory("results/slow5/{sample}")),
        completed=temp("results/slow5/{sample}.completed")
    log:
        "logs/slow5tools_fast5toslow5/{sample}.log"
    threads: config['threads']['slow5tools']
    # conda:
    #     "../envs/slow5tools.yaml"
    # the normal slow5tools will stuck by the uncompleted fast5 files
    container:
        "docker://btrspg/slow5tools:dev"
    shell:
        "slow5tools f2s {input.fast5} -d {output.outdir} -p {threads} --bad5 1 "
        "2>{log} && "
        "touch {output.completed}"

rule slow5tools_s2b:
    input:
        slow5="results/slow5/{sample}",
        completed="results/slow5/{sample}.completed"
    output:
        "results/blow5/{sample}.blow5"
    log:
        "logs/slow5tools_slow5toblow5/{sample}.log"
    threads: config['threads']['slow5tools']
    conda:
        "../envs/slow5tools.yaml"
    # container:
    #     "docker://btrspg/slow5tools:dev"
    shell:
        "slow5tools merge {input.slow5} -o {output} -t {threads} "
        "2>{log} "

