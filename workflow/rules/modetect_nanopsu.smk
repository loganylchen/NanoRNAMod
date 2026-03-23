rule nanopsu_predict:
    input:
        bam="{project}/results/alignments/{sample}_filtered.bam",
        bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        reference=config["reference"]["transcriptome_fasta"],
    output:
        outdir=KEEP_OR_NOT(directory("{project}/results/dataprep/{sample}_nanopsu_dataprep")),
    container:
        get_container("nanopsu")
    threads: get_threads("nanopsu", 4)
    resources:
        mem_mb = 1024 * 80
    priority: 10
    params:
        extra=config["params"].get("nanopsu", ""),
    log:
        "logs/{project}/nanopsu_predict/{sample}.log"
    benchmark:
        "benchmarks/{project}/{sample}.nanopsu.benchmark.txt"
    shell:
        "nanopsu predict "
        "--bam {input.bam} "
        "--reference {input.reference} "
        "--threads {threads} "
        "--output {output.outdir} "
        "{params.extra} "
        "2>{log}"
