rule uncompress_eventalign_full:
    input:
        completion="results/eventalign/{sample}_full.tsv.completed",
        eventalign="results/eventalign/{sample}_full.tsv.bz2",
    output:
        uc_eventalign = temp("results/eventalign/{sample}_full.tsv.bz2.tmp"),
        uc_completion = temp("results/eventalign/{sample}_full.tsv.completed.tmp")
    log:
        "logs/uncompress_eventalign_full/{sample}.log"
    benchmark:
        "benchmarks/{sample}.full_uncompress_eventalign.benchmark.txt"
    threads: 1
    shell:
        "bzip2 -dc {input.eventalign} > {output.uc_eventalign} && touch {output.uc_completion} 2>{log}"





rule nanocompore_collapse:
    input:
        eventalign="results/eventalign/{sample}_full.tsv.bz2.tmp",
        completion="results/eventalign/{sample}_full.tsv.completed.tmp"
    output:
        output=temp("results/nanocompore_eventalign_collapse/{sample}/{sample}_eventalign_collapse.tsv")
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
        dir=directory("results/nanocompore/{native}_{control}"),
        output_file="results/nanocompore/{native}_{control}/nanocompore_results.tsv"
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
        "--outpath {output.dir} "
        "--outprefix '' "
        "--overwrite  2>{log}"


