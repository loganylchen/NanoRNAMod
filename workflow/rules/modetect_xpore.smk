rule uncompress_eventalign:
    input:
        completion="{project}/results/eventalign/{sample}_simple.tsv.completed",
        eventalign="{project}/results/eventalign/{sample}_simple.tsv.bz2",
    output:
        uc_eventalign=temp("{project}/results/eventalign/{sample}_simple.tsv.bz2.tmp"),
        uc_completion=temp(
            "{project}/results/eventalign/{sample}_simple.tsv.completed.tmp"
        ),
    log:
        "logs/{project}/uncompress_eventalign/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.uncompress_eventalign.benchmark.txt"
    resources:
        mem_mb = 1024
    priority: 10
    threads: 1
    shell:
        "bzip2 -dc {input.eventalign} > {output.uc_eventalign} && touch {output.uc_completion} 2>{log}"


rule xpore_dataprep:
    input:
        completion="{project}/results/eventalign/{sample}_simple.tsv.completed.tmp",
        eventalign="{project}/results/eventalign/{sample}_simple.tsv.bz2.tmp",
        reference=config["reference"]["transcriptome_fasta"],
    output:
        KEEP_OR_NOT(directory("{project}/results/dataprep/{sample}_xpore_dataprep")),
    log:
        "logs/{project}/xpore_dataprep/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.xpore_dataprep.benchmark.txt"
    threads: config["threads"]["xpore"]
    priority: 10
    resources:
        mem_mb = 1024 * 50
    params:
        extra="",
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore dataprep "
        "--eventalign {input.eventalign} "
        "--transcript_fasta {input.reference} --n_processes {threads} "
        "{params.extra} "
        "--out_dir {output} 2>{log} "


rule xpore_config:
    input:
        control_dir="{project}/results/dataprep/{control}_xpore_dataprep",
        native_dir="{project}/results/dataprep/{native}_xpore_dataprep",
    output:
        conf="{project}/results/xpore/{native}_{control}.xpore_config.yaml",
    threads: 1
    resources:
        mem_mb = 1024
    params:
        "{project}/results/xpore/{native}_{control}",
    log:
        "logs/{project}/xpore_config/{native}_{control}.log",
    script:
        "../scripts/xpore_config.py"


rule xpore_run:
    input:
        "{project}/results/xpore/{native}_{control}.xpore_config.yaml",
        "{project}/results/dataprep/{control}_xpore_dataprep",
        "{project}/results/dataprep/{native}_xpore_dataprep",
    output:
        difftable=KEEP_OR_NOT(
            "{project}/results/xpore/{native}_{control}/diffmod.table"
        ),
    threads: config["threads"]["xpore"]
    resources:
        mem_mb = 1024 * 50
    priority: 10
    log:
        stdout="logs/{project}/xpore/{native}_{control}.log",
        stderr="logs/{project}/xpore/{native}_{control}.err",
    benchmark:
        "benchmarks/{project}/{native}_{control}.xpore.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore diffmod --config {input[0]} --n_processes {threads} 1>{log.stdout} 2>{log.stderr}"


rule xpore_postprocessing:
    input:
        "{project}/results/xpore/{native}_{control}/diffmod.table",
    output:
        KEEP_OR_NOT(
            "{project}/results/xpore/{native}_{control}/majority_direction_kmer_diffmod.table"
        ),
    threads: config["threads"]["xpore"]
    resources:
        mem_mb = 1024 * 50
    params:
        "{project}/results/xpore/{native}_{control}",
    threads: 1
    log:
        "logs/{project}/xpore_postprocessing/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.xpore_postprocessing.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore postprocessing --diffmod_dir {params}  2>{log}"
