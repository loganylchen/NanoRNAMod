rule directrm_config:
    input:
        bam="{project}/results/alignments/{sample}_filtered.bam",
        bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        reference=config["reference"]["transcriptome_fasta"]
    output:
        config="{project}/results/directrm/{sample}_directrm_config.yaml"
    container:
        get_container("directrm")
    threads: get_threads("directrm", 1)
    resources:
        mem_mb = 1024
    params:
        modifications="all",
        min_prob=0.7
    log:
        "logs/{project}/directrm_config/{sample}.log"
    script:
        "../scripts/directrm_config.py"


rule directrm_predict:
    input:
        bam="{project}/results/alignments/{sample}_filtered.bam",
        bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        reference=config["reference"]["transcriptome_fasta"],
        config="{project}/results/directrm/{sample}_directrm_config.yaml"
    output:
        predictions=temp("{project}/results/directrm/{sample}_predictions.tsv"),
        completion=touch("{project}/results/directrm/{sample}_predictions.completed")
    container:
        get_container("directrm")
    threads: get_threads("directrm", 4)
    resources:
        mem_mb = 1024 * 100
    priority: 10
    log:
        "logs/{project}/directrm_predict/{sample}.log"
    benchmark:
        "benchmarks/{project}/{sample}.directrm.benchmark.txt"
    params:
        extra=""
    shell:
        "directrm predict "
        "--bam {input.bam} "
        "--reference {input.reference} "
        "--config {input.config} "
        "--output {output.predictions} "
        "--threads {threads} "
        "{params.extra} "
        "2>{log} && touch {output.completion}"


rule directrm_postprocess:
    input:
        predictions="{project}/results/directrm/{sample}_predictions.completed",
        pred_file="{project}/results/directrm/{sample}_predictions.tsv"
    output:
        directory("{project}/results/dataprep/{sample}_directrm_dataprep")
    container:
        get_container("directrm")
    threads: get_threads("directrm", 1)
    resources:
        mem_mb = 1024 * 20
    priority: 10
    log:
        "logs/{project}/directrm_postprocess/{sample}.log"
    script:
        "../scripts/directrm_postprocess.py"
