rule xpore_dataprep:
    input:
        completion="results/eventalign/{sample}_xpore.tsv.completed",
        eventalign="results/eventalign/{sample}_xpore.tsv.gz",
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
        "gzip -dc {input.eventalign} > {input.eventalign}.xpore.tmp && xpore dataprep "
        "--eventalign {input.eventalign}.xpore.tmp "
        "--transcript_fasta {input.reference} "
        "--out_dir {output} 2>{log} && rm {input.eventalign}.xpore.tmp"


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

rule xpore_batch_config:
    input:
        control_dirs=expand("results/dataprep/{control}_xpore_dataprep",control=control_samples),
        native_dirs=expand("results/dataprep/{native}_xpore_dataprep",native=native_samples)
    output:
        conf="results/xpore/Batch_{native_list}_{control_list}.xpore_config.yaml"
    params:
        "results/xpore/Batch_{native_list}_{control_list}"
    log:
        "logs/xpore_config/{native_list}_{control_list}.log"
    script:
        "../scripts/batch_xpore_config.py"

rule xpore_run:
    input:
        "results/xpore/{native}_{control}.xpore_config.yaml"
    output:
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