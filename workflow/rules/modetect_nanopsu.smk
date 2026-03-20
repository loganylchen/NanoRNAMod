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
    conda:
        "../envs/nanopsu.yaml"
    shell:
        "nanopsu predict "
        "--bam {input.bam} "
        "--reference {input.reference} "
        "--threads {threads} "
        "--output {output.outdir} "
        "{params.extra} "
        "2>{log}"


rule nanopsu_postprocess:
    input:
        predictions="{project}/results/dataprep/{sample}_nanopsu_dataprep/predictions.tsv",
        dataprep_dir="{project}/results/dataprep/{sample}_nanopsu_dataprep",
    output:
        "{project}/results/modifications/nanopsu/{sample}/nanopsu_results.tsv",
    params:
        tool="nanopsu",
    container:
        get_container("default")
    threads: get_threads("default", 1)
    resources:
        mem_mb = 1024 * 20
    priority: 20
    log:
        "logs/{project}/nanopsu_postprocess/{sample}.log"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"
