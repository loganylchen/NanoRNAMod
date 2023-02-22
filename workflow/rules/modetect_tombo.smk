rule tombo_lsc_mod:
    input:
        completion1="results/multi_to_single/{native}/resquiggle_complete.txt",
        completion2="results/multi_to_single/{control}/resquiggle_complete.txt",
        control_dir="results/multi_to_single/{control}",
        native_dir="results/multi_to_single/{native}"
    output:
        "results/tombo_mod_lsc/{native}_{control}.lsc.tombo.stats"
    params:
        output_prefix="results/tombo_mod_lsc/{native}_{control}.lsc",
        extra=config['params']['tombo_lsc']
    log:
        "logs/tombo_mod_lsc/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.tombo_mod_lsc.benchmark.txt"
    threads: config['threads']['tombo']
    conda:
        "../envs/tombo.yaml"
    shell:
        "tombo detect_modifications level_sample_compare "
        "--fast5-basedirs {input.control_dir} "
        "--alternate-fast5-basedirs {input.native_dir} "
        "--statistics-file-basename {params.output_prefix} "
        "{params.extra} "
        "--processes {threads} 2>{log}"

rule tombo_mod_msc:
    input:
        completion1="results/multi_to_single/{native}/resquiggle_complete.txt",
        completion2="results/multi_to_single/{control}/resquiggle_complete.txt",
        control_dir="results/multi_to_single/{control}",
        native_dir="results/multi_to_single/{native}"
    output:
        transcript_res="results/tombo_mod_msc/{native}_{control}.msc.tombo.stats",
        molecule_res="results/tombo_mod_msc/{native}_{control}.msc.tombo.per_read_stats"
    params:
        output_prefix="results/tombo_mod_msc/{native}_{control}.msc",
        extra=config['params']['tombo_msc']
    log:
        "logs/tombo_mod_msc/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.tombo_mod_msc.benchmark.txt"
    threads: config['threads']['tombo']
    conda:
        "../envs/tombo.yaml"

    shell:
        "tombo detect_modifications model_sample_compare "
        "--fast5-basedirs {input.native_dir} "
        "--control-fast5-basedirs {input.control_dir} "
        "{params.extra} "
        "--per-read-statistics-basename {params.output_prefix} "
        "--processes {threads} "
        "--statistics-file-basename {params.output_prefix}  2>{log}"

rule tombo_resquiggle:
    input:
        fast5="results/multi_to_single/{sample}",
        reference=config['reference']['transcriptome_fasta'],
        completion="results/multi_to_single/{sample}/m2s_complete.txt"
    output:
        completion="results/multi_to_single/{sample}/resquiggle_complete.txt"
    log:
        "logs/tombo_resquiggle/{sample}.log"
    benchmark:
        "benchmarks/{sample}.tombo_resquiggle.benchmark.txt"
    params:
        extra=config['params']['tombo_resquiggle']
    threads: config['threads']['tombo']
    conda:
        "../envs/tombo.yaml"
    shell:
        "tombo resquiggle "
        "{input.fast5} "
        "{input.reference} "
        "--processes {threads} "
        "{params.extra}  2>{log} && echo `date` > {output.completion} "