rule rnano_config:
    input:
        bam="{project}/results/alignments/{sample}_filtered.bam",
        bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        reference=config["reference"]["transcriptome_fasta"]
    output:
        config="{project}/results/rnano/{sample}_rnano_config.yaml"
    container:
        get_container("rnano")
    threads: get_threads("rnano", 1)
    resources:
        mem_mb = 1024
    params:
        model="pretrained",
        gpu=False
    log:
        "logs/{project}/rnano_config/{sample}.log"
    script:
        "../scripts/rnano_config.py"


rule rnano_predict:
    input:
        bam="{project}/results/alignments/{sample}_filtered.bam",
        bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        reference=config["reference"]["transcriptome_fasta"],
        config="{project}/results/rnano/{sample}_rnano_config.yaml"
    output:
        predictions=temp("{project}/results/rnano/{sample}_predictions.tsv"),
        completion=touch("{project}/results/rnano/{sample}_predictions.completed")
    container:
        get_container("rnano")
    threads: get_threads("rnano", 4)
    resources:
        mem_mb = 1024 * 100
    priority: 10
    log:
        "logs/{project}/rnano_predict/{sample}.log"
    benchmark:
        "benchmarks/{project}/{sample}.rnano.benchmark.txt"
    params:
        extra=""
    shell:
        "rnano predict "
        "--bam {input.bam} "
        "--reference {input.reference} "
        "--config {input.config} "
        "--threads {threads} "
        "{params.extra} "
        "--output {output.predictions} "
        "2>{log} && touch {output.completion}"


rule rnano_postprocess:
    input:
        predictions="{project}/results/rnano/{sample}_predictions.completed",
        pred_file="{project}/results/rnano/{sample}_predictions.tsv"
    output:
        directory("{project}/results/dataprep/{sample}_rnano_dataprep")
    container:
        get_container("rnano")
    threads: get_threads("rnano", 1)
    resources:
        mem_mb = 1024 * 20
    priority: 10
    log:
        "logs/{project}/rnano_postprocess/{sample}.log"
    script:
        "../scripts/rnano_postprocess.py"
