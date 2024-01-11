rule xpore_dataprep:
    input:
        completion="results/eventalign/{sample}_xpore.tsv.completed",
        eventalign="results/eventalign/{sample}_xpore.tsv.gz",
        reference=config['reference']['transcriptome_fasta'],
    output:
        directory("results/dataprep/{sample}_xpore_dataprep")
    log:
        "logs/xpore_dataprep/{sample}.log"
    benchmark:
        "benchmarks/{sample}.xpore_dataprep.benchmark.txt"
    params:
        extra=''
    conda:
        "../envs/xpore.yaml"
    shell:
        "gzip -dc {input.eventalign} > {input.eventalign}.xpore.tmp && xpore dataprep "
        "--eventalign {input.eventalign}.xpore.tmp "
        "--transcript_fasta {input.reference} "
        "{params.extra} "
        "--out_dir {output} 2>{log} && rm {input.eventalign}.xpore.tmp"



rule xpore_dataprep_genome:
    input:
        completion="results/eventalign/{sample}_xpore.tsv.completed",
        eventalign="results/eventalign/{sample}_xpore.tsv.gz",
        transcriptome=config['reference']['transcriptome_fasta'],
        gtf=config['reference']['transcriptome_gtf']
    output:
        directory("results/dataprep/{sample}_xpore_dataprep_gnome")
    log:
        "logs/xpore_dataprep/{sample}_genome.log"
    benchmark:
        "benchmarks/{sample}.xpore_dataprep_genome.benchmark.txt"
    params:
        extra=''
    conda:
        "../envs/xpore.yaml"
    shell:
        "gzip -dc {input.eventalign} > {input.eventalign}.xpore.tmp && xpore dataprep "
        "--eventalign {input.eventalign}.xpore.tmp "
        "--transcript_fasta {input.transcriptome} "
        "--gtf_or_gff {input.gtf} "
        "--genome "
        "{params.extra} "
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

rule xpore_config_genome:
    input:
        control_dir="results/dataprep/{control}_xpore_dataprep_genome",
        native_dir="results/dataprep/{native}_xpore_dataprep_genome"
    output:
        conf="results/xpore/{native}_{control}.xpore_config_genome.yaml"
    params:
        "results/xpore/{native}_{control}_genome"
    log:
        "logs/xpore_config/{native}_{control}_genome.log"
    script:
        "../scripts/xpore_config.py"

rule xpore_group_config:
    input:
        control_dirs=expand("results/dataprep/{control}_xpore_dataprep",control=control_samples),
        native_dirs=expand("results/dataprep/{native}_xpore_dataprep",native=native_samples)
    output:
        conf="results/xpore/Groups/{native_list}_{control_list}.xpore_config.yaml"
    params:
        "results/xpore/Groups/{native_list}_{control_list}"
    log:
        "logs/xpore_config/Groups_{native_list}_{control_list}.log"
    script:
        "../scripts/batch_xpore_config.py"

rule xpore_group_config_genome:
    input:
        control_dirs=expand("results/dataprep/{control}_xpore_dataprep_genome",control=control_samples),
        native_dirs=expand("results/dataprep/{native}_xpore_dataprep_genome",native=native_samples)
    output:
        conf="results/xpore/Groups/{native_list}_{control_list}.xpore_config_genome.yaml"
    params:
        "results/xpore/Groups/{native_list}_{control_list}_genome"
    log:
        "logs/xpore_config/Groups_{native_list}_{control_list}_genome.log"
    script:
        "../scripts/batch_xpore_config.py"

rule xpore_group_run:
    input:
        "results/xpore/Groups/{native_list}_{control_list}.xpore_config.yaml"
    output:
        difftable="results/xpore/Groups/{native_list}_{control_list}/diffmod.table"
    log:
        "logs/xpore/Groups_{native_list}_{control_list}.log"
    benchmark:
        "benchmarks/Groups_{native_list}_{control_list}.xpore.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore diffmod --config {input} 2>{log}"


rule xpore_group_run_genome:
    input:
        "results/xpore/Groups/{native_list}_{control_list}.xpore_config_genome.yaml"
    output:
        difftable="results/xpore/Groups/{native_list}_{control_list}_genome/diffmod.table"
    log:
        "logs/xpore/Groups_{native_list}_{control_list}_genome.log"
    benchmark:
        "benchmarks/Groups_{native_list}_{control_list}.xpore_genome.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore diffmod --config {input} 2>{log}"

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

rule xpore_run_genome:
    input:
        "results/xpore/{native}_{control}.xpore_config_genome.yaml"
    output:
        difftable="results/xpore/{native}_{control}_genome/diffmod.table"
    log:
        "logs/xpore/{native}_{control}_genome.log"
    benchmark:
        "benchmarks/{native}_{control}.xpore_genome.benchmark.txt"
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

rule xpore_postprocessing_genome:
    input:
        "results/xpore/{native}_{control}_genome/diffmod.table"
    output:
        "results/xpore/{native}_{control}_genome/majority_direction_kmer_diffmod.table"
    params:
        "results/xpore/{native}_{control}_genome"
    log:
        "logs/xpore_postprocessing/{native}_{control}_genome.log"
    benchmark:
        "benchmarks/{native}_{control}.xpore_postprocessing_genome.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore postprocessing --diffmod_dir {params}  2>{log}"

rule xpore_group_postprocessing:
    input:
        "results/xpore/Groups/{native_list}_{control_list}/diffmod.table"
    output:
        "results/xpore/Groups/{native_list}_{control_list}/majority_direction_kmer_diffmod.table"
    params:
        "results/xpore/Groups/{native_list}_{control_list}"
    log:
        "logs/xpore_postprocessing/Groups_{native_list}_{control_list}.log"
    benchmark:
        "benchmarks/Groups_{native_list}_{control_list}.xpore_postprocessing.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore postprocessing --diffmod_dir {params}  2>{log}"


rule xpore_group_postprocessing_genome:
    input:
        "results/xpore/Groups/{native_list}_{control_list}_genome/diffmod.table"
    output:
        "results/xpore/Groups/{native_list}_{control_list}_genome/majority_direction_kmer_diffmod.table"
    params:
        "results/xpore/Groups/{native_list}_{control_list}_genome"
    log:
        "logs/xpore_postprocessing/Groups_{native_list}_{control_list}_genome.log"
    benchmark:
        "benchmarks/Groups_{native_list}_{control_list}.xpore_postprocessing_genome.benchmark.txt"
    conda:
        "../envs/xpore.yaml"
    shell:
        "xpore postprocessing --diffmod_dir {params}  2>{log}"