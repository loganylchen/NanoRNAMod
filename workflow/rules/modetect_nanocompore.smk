rule nanocompore_collapse:
    input:
        eventalign="results/eventalign/{sample}_nanocompore.tsv.gz",
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

rule nanocompore_collapse_sampled:
    input:
        eventalign="results/eventalign/{sample}_nanocompore_{sample_size}.tsv.gz",
        completion="results/eventalign/{sample}_nanocompore_{sample_size}.tsv.completed"
    output:
        output="results/nanocompore_eventalign_collapse/{sample}_{sample_size}/{sample}_{sample_size}_eventalign_collapse.tsv"
    params:
        prefix="{sample}_{sample_size}",
        dir="results/nanocompore_eventalign_collapse/{sample}_{sample_size}",
    log:
        "logs/nanocompore_collapse/{sample}_{sample_size}.log"
    benchmark:
        "benchmarks/{sample}.nanocompore_collapse_{sample_size}.benchmark.txt"
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

rule nanocompore_sampled:
    input:
        control_file="results/nanocompore_eventalign_collapse/{control}_{sample_size}/{control}_{sample_size}_eventalign_collapse.tsv",
        native_file="results/nanocompore_eventalign_collapse/{native}_{sample_size}/{native}_{sample_size}_eventalign_collapse.tsv",
        reference=config['reference']['transcriptome_fasta']
    output:
        directory("results/nanocompore/{native}_{control}-{sample_size}"),
        "results/nanocompore/{native}_{control}-{sample_size}/{native}_{control}_{sample_size}nanocompore_results.tsv"
    params:
        prefix="{native}_{control}_{sample_size}",
        extra=config['params']['nanocompore']
    log:
        stdout="logs/nanocompore/{native}_{control}_{sample_size}.log"
    benchmark:
        "benchmarks/{native}_{control}.nanocompore_{sample_size}.benchmark.txt"
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

rule nanocompore_group:
    input:
        control_file=expand("results/nanocompore_eventalign_collapse/{sample}/{sample}_eventalign_collapse.tsv",sample=control_samples),
        native_file=expand("results/nanocompore_eventalign_collapse/{sample}/{sample}_eventalign_collapse.tsv",sample=native_samples),
        reference=config['reference']['transcriptome_fasta']
    output:
        directory("results/nanocompore/Group_{native_list}_{control_list}")
    params:
        prefix="Group_{native_list}_{control_list}",
        extra=config['params']['nanocompore'],
        control_files=get_nanocompore_list(control_samples),
        native_files=get_nanocompore_list(native_samples),
    log:
        stdout="logs/nanocompore/Group_{native_list}_{control_list}.log"
    benchmark:
        "benchmarks/Group_{native_list}_{control_list}.nanocompore.benchmark.txt"
    conda:
        "../envs/nanocompore.yaml"
    shell:
        "nanocompore sampcomp "
        "--file_list1 {params.control_files} "
        "--file_list2 {params.native_files} "
        "--label1 Control "
        "--label2 Native "
        "{params.extra} "
        "--fasta {input.reference} "
        "--outpath {output} "
        "--outprefix {params.prefix} "
        "--overwrite  2>{log}"