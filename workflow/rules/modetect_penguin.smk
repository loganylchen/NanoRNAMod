rule penguin_predict:
    input:
        bam="{project}/results/alignments/{sample}_filtered.bam",
        bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        reference=config["reference"]["transcriptome_fasta"],
    output:
        outdir=KEEP_OR_NOT(directory("{project}/results/dataprep/{sample}_penguin_dataprep")),
    container:
        get_container("penguin")
    threads: get_threads("penguin", 4)
    resources:
        mem_mb = 1024 * 50
    priority: 10
    params:
        extra=config["params"].get("penguin", ""),
    log:
        "logs/{project}/penguin_predict/{sample}.log"
    benchmark:
        "benchmarks/{project}/{sample}.penguin.benchmark.txt"
    shell:
        "penguin predict "
        "--bam {input.bam} "
        "--reference {input.reference} "
        "--threads {threads} "
        "--output {output.outdir} "
        "{params.extra} "
        "2>{log}"
