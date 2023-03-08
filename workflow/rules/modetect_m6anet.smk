rule dataprep_m6anet:
    input:
        eventalign="results/eventalign/{sample}_xpore.tsv.gz",
        completion="results/eventalign/{sample}_xpore.tsv.completed"
    output:
        directory("results/dataprep/{sample}_m6anet_dataprep")
    log:
        "logs/m6anet_dataprep/{sample}.log"
    benchmark:
        "benchmarks/{sample}.m6anet_dataprep.benchmark.txt"
    threads: config['threads']['m6anet']
    conda:
        "../envs/m6anet.yaml"
    shell:
        "gzip -dc {input.eventalign} > {input.eventalign}.m6anet.tmp &&"
        "m6anet-dataprep "
        "--eventalign {input.eventalign}.m6anet.tmp "
        "--n_processes {threads} "
        "--out_dir {output} 2>{log} && rm {input.eventalign}.m6anet.tmp"

rule m6anet_inference:
    input:
        "results/dataprep/{sample}_m6anet_dataprep"
    output:
        dir=directory("results/m6anet/{sample}"),
        result="results/m6anet/{sample}/data.result.csv.gz"
    log:
        "logs/m6anet_inference/{sample}.log"
    benchmark:
        "benchmarks/{sample}.m6anet_inference.benchmark.txt"
    threads: config['threads']['m6anet']
    conda:
        "../envs/m6anet.yaml"
    shell:
        "m6anet-run_inference "
        "--input_dir {input} "
        "--out_dir {output.dir} "
        "--infer_mod_rate "
        "--n_processes {threads}  2>{log}"
