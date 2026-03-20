rule tandemmod_config:
    input:
        bam="{project}/results/alignments/{sample}_filtered.bam",
        bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        reference=config["reference"]["transcriptome_fasta"]
    output:
        config="{project}/results/tandemmod/{sample}_tandemmod_config.yaml"
    container:
        get_container("tandemmod")
    threads: get_threads("tandemmod", 1)
    resources:
        mem_mb = 1024
    params:
        model_type="multi",
        threshold=0.5
    log:
        "logs/{project}/tandemmod_config/{sample}.log"
    script:
        "../scripts/tandemmod_config.py"


rule tandemmod_predict:
    input:
        bam="{project}/results/alignments/{sample}_filtered.bam",
        bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        reference=config["reference"]["transcriptome_fasta"],
        config="{project}/results/tandemmod/{sample}_tandemmod_config.yaml"
    output:
        predictions=temp("{project}/results/tandemmod/{sample}_predictions.tsv"),
        completion=touch("{project}/results/tandemmod/{sample}_predictions.completed")
    container:
        get_container("tandemmod")
    threads: get_threads("tandemmod", 4)
    resources:
        mem_mb = 1024 * 100
    priority: 10
    log:
        "logs/{project}/tandemmod_predict/{sample}.log"
    benchmark:
        "benchmarks/{project}/{sample}.tandemmod.benchmark.txt"
    params:
        extra=""
    conda:
        "../envs/tandemmod.yaml"
    shell:
        "tandemmod predict "
        "--bam {input.bam} "
        "--reference {input.reference} "
        "--config {input.config} "
        "--threads {threads} "
        "{params.extra} "
        "--output {output.predictions} "
        "2>{log} && touch {output.completion}"


rule tandemmod_postprocess:
    input:
        predictions="{project}/results/tandemmod/{sample}_predictions.completed",
        pred_file="{project}/results/tandemmod/{sample}_predictions.tsv"
    output:
        directory("{project}/results/dataprep/{sample}_tandemmod_dataprep")
    container:
        get_container("tandemmod")
    threads: get_threads("tandemmod", 1)
    resources:
        mem_mb = 1024 * 20
    priority: 10
    log:
        "logs/{project}/tandemmod_postprocess/{sample}.log"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/tandemmod_postprocess.py"
