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
    container:
        get_container("default")
    resources:
        mem_mb = 1024 * 10
    priority: 90
    threads: get_threads("default", 1)
    shell:
        "bzip2 -dc {input.eventalign} > {output.uc_eventalign} && touch {output.uc_completion} 2>{log}"


rule xpore_dataprep:
    input:
        completion="{project}/results/eventalign/{sample}_simple.tsv.completed.tmp",
        eventalign="{project}/results/eventalign/{sample}_simple.tsv.bz2.tmp",
        reference=config["reference"]["transcriptome_fasta"],
    output:
        temp(directory("{project}/results/dataprep/{sample}_xpore_dataprep")),
    log:
        "logs/{project}/xpore_dataprep/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.xpore_dataprep.benchmark.txt"
    container:
        get_container("xpore")
    threads: get_threads("xpore", 4)
    priority: 91
    resources:
        mem_mb = 1024 * 50
    params:
        extra="",
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
        conf=temp("{project}/results/xpore/{native}_{control}.xpore_config.yaml"),
    container:
        get_container("xpore")
    threads: get_threads("default", 1)
    resources:
        mem_mb = 1024 * 10
    priority: 92
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
        difftable=temp(
            "{project}/results/xpore/{native}_{control}/diffmod.table"
        ),
    container:
        get_container("xpore")
    threads: get_threads("xpore", 4)
    resources:
        mem_mb = 1024 * 50
    priority: 93
    log:
        "logs/{project}/xpore/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.xpore.benchmark.txt"
    shell:
        "xpore diffmod --config {input[0]} --n_processes {threads} 1>{log} 2>&1; "
        "status=$?; "
        "if [ $status -eq 0 ] && [ ! -f {output.difftable} ]; then "
        "  mkdir -p $(dirname {output.difftable}) && touch {output.difftable}; "
        "fi; "
        "exit $status"


rule xpore_postprocessing:
    input:
        "{project}/results/xpore/{native}_{control}/diffmod.table",
    output:
        temp(
            "{project}/results/xpore/{native}_{control}/majority_direction_kmer_diffmod.table"
        ),
    container:
        get_container("xpore")
    threads: get_threads("default", 1)
    resources:
        mem_mb = 1024 * 50
    priority: 94
    params:
        "{project}/results/xpore/{native}_{control}",
    log:
        "logs/{project}/xpore_postprocessing/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.xpore_postprocessing.benchmark.txt"
    shell:
        "xpore postprocessing --diffmod_dir {params}  2>{log}"
