rule nanocompore_collapse:
    input:
        eventalign="results/eventalign/{sample}_nanocompore.tsv",
        completion="results/eventalign/{sample}_nanocompore.tsv.completed"
    output:
        output="results/nanocompore_eventalign_collapse/{sample}/{sample}_eventalign_collapse.tsv"
    params:
        prefix="{sample}",
        dir="results/nanocompore_eventalign_collapse/{sample}",
    log:
        "logs/nanocompore_collapse/{sample}.log"
    benchmark:
        "benchmarks/{sample}.nanocompore_collapse.benchmark.txt"
    conda:
        "../envs/nanocompore.yaml"
    threads: config['threads']['nanocompore']
    shell:
        "nanocompore eventalign_collapse "
        "-i {input.eventalign} "
        "--outpath {params.dir} "
        "--outprefix {params.prefix} "
        "--overwrite "
        "--nthreads {threads} 2>{log}"

rule nanocompore:
    input:
        control_file="results/nanocompore_eventalign_collapse/{control}/{control}_eventalign_collapse.tsv",
        native_file="results/nanocompore_eventalign_collapse/{native}/{native}_eventalign_collapse.tsv",
        reference=config['reference']['transcriptome_fasta']
    output:
        directory("results/nanocompore/{native}_{control}")
    params:
        prefix="{native}_{control}",
        extra=config['params']['nanocompore']
    log:
        stdout="logs/nanocompore/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.nanocompore.benchmark.txt"
    conda:
        "../envs/nanocompore.yaml"
    shell:
        "nanocompore sampcomp "
        "--file_list1 {input.control_file} "
        "--file_list2 {input.native_file} "
        "--label1 Control "
        "--label2 Native "
        "{params.extra} "
        "--fasta {input.reference} "
        "--outpath {output} "
        "--outprefix {params.prefix} "
        "--overwrite  2>{log}"