rule nanomud_predict:
    input:
        bam="{project}/results/alignments/{sample}_filtered.bam",
        bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        reference=config["reference"]["transcriptome_fasta"],
    output:
        outdir=KEEP_OR_NOT(directory("{project}/results/dataprep/{sample}_nanomud_dataprep")),
    container:
        get_container("nanomud")
    threads: get_threads("nanomud", 4)
    resources:
        mem_mb = 1024 * 80
    priority: 10
    params:
        mods=config["params"].get("nanomud_mods", "Psi,m1Psi"),
        extra=config["params"].get("nanomud", ""),
    log:
        "logs/{project}/nanomud_predict/{sample}.log"
    benchmark:
        "benchmarks/{project}/{sample}.nanomud.benchmark.txt"
    conda:
        "../envs/nanomud.yaml"
    shell:
        "nanomud predict "
        "--bam {input.bam} "
        "--reference {input.reference} "
        "--mods {params.mods} "
        "--threads {threads} "
        "--output {output.outdir} "
        "{params.extra} "
        "2>{log}"
