rule dataprep_m6anet:
    input:
        completion="results/eventalign/{sample}_xpore.tsv.completed.tmp",
        eventalign="results/eventalign/{sample}_xpore.tsv.gz.tmp",
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
        ""
        "m6anet dataprep "
        "--eventalign {input.eventalign} "
        "--n_processes {threads} --compress "
        "--out_dir {output} 2>{log} "

rule m6anet_inference:
    input:
        "results/dataprep/{sample}_m6anet_dataprep"
    output:
        dir=directory("results/m6anet/{sample}"),
        result="results/m6anet/{sample}/data.site_proba.csv"
    log:
        "logs/m6anet_inference/{sample}.log"
    benchmark:
        "benchmarks/{sample}.m6anet_inference.benchmark.txt"
    threads: config['threads']['m6anet']
    conda:
        "../envs/m6anet.yaml"
    shell:
        "m6anet inference "
        "--input_dir {input} "
        "--out_dir {output.dir} "
        "--n_processes {threads}  2>{log}"
