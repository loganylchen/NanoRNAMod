rule xpore_dataprep:
    input:
        completion="results/eventalign/{sample}_xpore.tsv.completed",
        eventalign="results/eventalign/{sample}_xpore.tsv",
        reference=config['reference']['transcriptome_fasta']
    output:
        directory("results/dataprep/{sample}_xpore_dataprep")
    log:
        "logs/xpore_dataprep/{sample}.log"
    benchmark:
        "benchmarks/{sample}.xpore_dataprep.benchmark.txt"

    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore dataprep "
        "--eventalign {input.eventalign} "
        "--transcript_fasta {input.reference} "
        "--out_dir {output} 2>{log}"


rule xpore_config:
    input:
        control_dir="results/dataprep/{control}_xpore_dataprep",
        native_dir="results/dataprep/{native}_xpore_dataprep"
    output:
        conf="results/xpore/{native}_{control}.xpore_config.yaml"
    params:
        "results/xpore/{native}_{control}"
    log:
        "logs/xpore_config/{native}_{control}.log"
    script:
        "../scripts/xpore_config.py"

rule xpore_run:
    input:
        "results/xpore/{native}_{control}.xpore_config.yaml"
    output:
        res_dir="results/xpore/{native}_{control}",
        difftable="results/xpore/{native}_{control}/diffmod.table"
    log:
        "logs/xpore/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.xpore.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore diffmod --config {input} 2>{log}"

rule xpore_postprocessing:
    input:
        "results/xpore/{native}_{control}/diffmod.table"
    output:
        "results/xpore/{native}_{control}/majority_direction_kmer_diffmod.table"
    params:
        "results/xpore/{native}_{control}"
    log:
        "logs/xpore_postprocessing/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.xpore_postprocessing.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore postprocessing --diffmod_dir {params}  2>{log}"