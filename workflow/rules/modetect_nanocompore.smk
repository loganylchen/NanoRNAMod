rule uncompress_eventalign_full:
    input:
        completion="{project}/results/eventalign/{sample}_full.tsv.completed",
        eventalign="{project}/results/eventalign/{sample}_full.tsv.bz2",
    output:
        uc_eventalign=temp("{project}/results/eventalign/{sample}_full.tsv.bz2.tmp"),
        uc_completion=temp(
            "{project}/results/eventalign/{sample}_full.tsv.completed.tmp"
        ),
    log:
        "logs/{project}/uncompress_eventalign_full/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.full_uncompress_eventalign.benchmark.txt"
    threads: 1
    shell:
        "bzip2 -dc {input.eventalign} > {output.uc_eventalign} && touch {output.uc_completion} 2>{log}"


rule nanocompore_collapse:
    input:
        eventalign="{project}/results/eventalign/{sample}_full.tsv.bz2.tmp",
        completion="{project}/results/eventalign/{sample}_full.tsv.completed.tmp",
    output:
        output=temp(
            "{project}/results/nanocompore_eventalign_collapse/{sample}/{sample}_eventalign_collapse.tsv"
        ),
    params:
        prefix="{sample}",
        dir="{project}/results/nanocompore_eventalign_collapse/{sample}",
    log:
        "logs/{project}/nanocompore_collapse/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.nanocompore_collapse.benchmark.txt"
    conda:
        "../envs/nanocompore.yaml"
    threads: config["threads"]["nanocompore"]
    shell:
        "nanocompore eventalign_collapse "
        "-i {input.eventalign} "
        "--outpath {params.dir} "
        "--outprefix {params.prefix} "
        "--overwrite "
        "--nthreads {threads} 2>{log}"


rule nanocompore:
    input:
        control_file="{project}/results/nanocompore_eventalign_collapse/{control}/{control}_eventalign_collapse.tsv",
        native_file="{project}/results/nanocompore_eventalign_collapse/{native}/{native}_eventalign_collapse.tsv",
        reference=config["reference"]["transcriptome_fasta"],
    output:
        dir=directory("{project}/results/nanocompore/{native}_{control}"),
        output_file="{project}/results/nanocompore/{native}_{control}/nanocompore_results.tsv",
    params:
        prefix="{native}_{control}",
        extra=config["params"]["nanocompore"],
    log:
        stdout="logs/{project}/nanocompore/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.nanocompore.benchmark.txt"
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
