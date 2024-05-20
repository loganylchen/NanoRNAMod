rule uncompress_eventalign:
    input:
        completion="results/eventalign/{sample}_simple.tsv.completed",
        eventalign="results/eventalign/{sample}_simple.tsv.bz2",
    output:
        uc_eventalign = temp("results/eventalign/{sample}_simple.tsv.bz2.tmp"),
        uc_completion = temp("results/eventalign/{sample}_simple.tsv.completed.tmp")
    log:
        "logs/uncompress_eventalign/{sample}.log"
    benchmark:
        "benchmarks/{sample}.uncompress_eventalign.benchmark.txt"
    threads: 1
    shell:
        "bzip2 -dc {input.eventalign} > {output.uc_eventalign} && touch {output.uc_completion} 2>{log}"


rule uncompress_eventalign_sampled:
    input:
        completion="results/eventalign/{sample}_simple_{sample_size}_{n}.tsv.completed",
        eventalign="results/eventalign/{sample}_simple_{sample_size}_{n}.tsv.bz2",
    output:
        uc_eventalign = temp("results/eventalign/{sample}_simple_{sample_size}_{n}.tsv.bz2.tmp"),
        uc_completion = temp("results/eventalign/{sample}_simple_{sample_size}_{n}.tsv.completed.tmp")
    log:
        "logs/uncompress_eventalign_{sample_size}_{n}/{sample}.log"
    benchmark:
        "benchmarks/{sample}.uncompress_eventalign_{sample_size}_{n}.benchmark.txt"
    threads: 1
    shell:
        "gzip -dc {input.eventalign} > {output.uc_eventalign} && touch {output.uc_completion} 2>{log}"

rule xpore_dataprep:
    input:
        completion="results/eventalign/{sample}_simple.tsv.completed.tmp",
        eventalign="results/eventalign/{sample}_simple.tsv.bz2.tmp",
        reference=config['reference']['transcriptome_fasta'],
    output:
        directory("results/dataprep/{sample}_xpore_dataprep")
    log:
        "logs/xpore_dataprep/{sample}.log"
    benchmark:
        "benchmarks/{sample}.xpore_dataprep.benchmark.txt"
    threads: config['threads']['xpore']
    params:
        extra=''
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore dataprep "
        "--eventalign {input.eventalign} "
        "--transcript_fasta {input.reference} --n_processes {threads} "
        "{params.extra} "
        "--out_dir {output} 2>{log} "


rule xpore_dataprep_sampled:
    input:
        completion="results/eventalign/{sample}_simple_{sample_size}_{n}.tsv.completed.tmp",
        eventalign="results/eventalign/{sample}_simple_{sample_size}_{n}.tsv.bz2.tmp",
        reference=config['reference']['transcriptome_fasta'],
    output:
        directory("results/dataprep/{sample}_xpore_dataprep_{sample_size}_{n}")
    log:
        "logs/xpore_dataprep_{sample_size}_{n}/{sample}.log"
    benchmark:
        "benchmarks/{sample}.xpore_dataprep_{sample_size}_{n}.benchmark.txt"
    threads: config['threads']['xpore']
    params:
        extra=''
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
        control_dir="results/dataprep/{control}_xpore_dataprep",
        native_dir="results/dataprep/{native}_xpore_dataprep"
    output:
        conf="results/xpore/{native}_{control}.xpore_config.yaml"
    threads: 1
    params:
        "results/xpore/{native}_{control}"
    log:
        "logs/xpore_config/{native}_{control}.log"
    script:
        "../scripts/xpore_config.py"


rule xpore_config_sampled:
    input:
        control_dir="results/dataprep/{control}_xpore_dataprep_{sample_size}_{n}",
        native_dir="results/dataprep/{native}_xpore_dataprep_{sample_size}_{n}"
    output:
        conf="results/xpore/{native}_{control}-{sample_size}-{n}.xpore_config.yaml"
    threads: 1
    params:
        "results/xpore/{native}_{control}-{sample_size}-{n}"
    log:
        "logs/xpore_config/{native}_{control}_{sample_size}_{n}.log"
    script:
        "../scripts/xpore_config.py"





rule xpore_run:
    input:
        "results/xpore/{native}_{control}.xpore_config.yaml"
    output:
        difftable="results/xpore/{native}_{control}/diffmod.table"
    threads: config['threads']['xpore']
    log:
        "logs/xpore/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.xpore.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore diffmod --config {input} --n_processes {threads} 2>{log}"

rule xpore_run_sampled:
    input:
        "results/xpore/{native}_{control}-{sample_size}-{n}.xpore_config.yaml"
    output:
        difftable="results/xpore/{native}_{control}-{sample_size}-{n}/diffmod.table"
    threads: config['threads']['xpore']
    log:
        "logs/xpore/{native}_{control}_{sample_size}_{n}.log"
    benchmark:
        "benchmarks/{native}_{control}_{sample_size}_{n}.xpore.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore diffmod --config {input} --n_processes {threads} 2>{log}"



rule xpore_postprocessing:
    input:
        "results/xpore/{native}_{control}/diffmod.table"
    output:
        "results/xpore/{native}_{control}/majority_direction_kmer_diffmod.table"
    threads: config['threads']['xpore']
    params:
        "results/xpore/{native}_{control}"
    threads: 1
    log:
        "logs/xpore_postprocessing/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.xpore_postprocessing.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore postprocessing --diffmod_dir {params}  2>{log}"

rule xpore_postprocessing_sampled:
    input:
        "results/xpore/{native}_{control}-{sample_size}-{n}/diffmod.table"
    output:
        "results/xpore/{native}_{control}-{sample_size}-{n}/majority_direction_kmer_diffmod.table"
    threads: config['threads']['xpore']
    params:
        "results/xpore/{native}_{control}-{sample_size}-{n}"
    threads: 1
    log:
        "logs/xpore_postprocessing/{native}_{control}_{sample_size}_{n}.log"
    benchmark:
        "benchmarks/{native}_{control}_{sample_size}_{n}.xpore_postprocessing.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore postprocessing --diffmod_dir {params}  2>{log}"


