
rule qc_nanoplot:
    input:
        fastq="{project}/results/fastq/{sample}.fq.gz",
    output:
        report="{project}/results/qc/{sample}/NanoStats.txt",
    params:
        output_dir=directory("{project}/results/qc/{sample}"),
    benchmark:
        "benchmarks/{project}/{sample}.nanoplot.benchmark.txt"
    log:
        "logs/{project}/qc/{sample}_nanoplot.log",
    conda:
        "../envs/nanoplot.yaml"
    resources:
        mem_mb=1024 * 20,
    priority: 50
    threads: config["threads"]["nanoplot"]
    shell:
        "NanoPlot --raw "
        "--outdir {params.output_dir} "
        "--tsv_stats --info_in_report "
        "-t {threads} --fastq {input.fastq} > {log}"
