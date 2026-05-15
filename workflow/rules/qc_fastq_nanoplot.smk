
rule qc_nanoplot:
    input:
        fastq="{project}/results/fastq/{sample}.fq.gz",
    output:
        report="{project}/results/qc/{sample}/NanoStats.txt",
    params:
        output_dir=directory("{project}/results/qc/{sample}"),
    container:
        get_container("nanoplot")
    benchmark:
        "benchmarks/{project}/{sample}.nanoplot.benchmark.txt"
    log:
        "logs/{project}/qc/{sample}_nanoplot.log",
    resources:
        mem_mb=1024 * 20,
    priority: 50
    threads: get_threads("nanoplot", 4)
    shell:
        "NanoPlot --raw "
        "--outdir {params.output_dir} "
        "--tsv_stats --info_in_report "
        "-t {threads} --fastq {input.fastq} > {log}"
