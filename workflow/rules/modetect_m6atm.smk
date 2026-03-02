rule m6atm_config:
    input:
        bam="{project}/results/bam/{sample}.bam",
        reference=config["reference"]["transcriptome_fasta"]
    output:
        config="{project}/results/m6atm/{sample}_m6atm_config.yaml"
    threads: 1
    resources:
        mem_mb = 1024
    params:
        threshold=0.6
        stoichiometry=True
    log:
        "logs/{project}/m6atm_config/{sample}.log"
    script:
        "../scripts/m6atm_config.py"


rule m6atm_predict:
    input:
        bam="{project}/results/bam/{sample}.bam",
        reference=config["reference"]["transcriptome_fasta"],
        config="{project}/results/m6atm/{sample}_m6atm_config.yaml"
    output:
        predictions=temp("{project}/results/m6atm/{sample}_predictions.tsv"),
        completion=touch("{project}/results/m6atm/{sample}_predictions.completed")
    threads: config["threads"]["m6atm"]
    resources:
        mem_mb = 1024 * 80
    priority: 10
    log:
        "logs/{project}/m6atm_predict/{sample}.log"
    benchmark:
        "benchmarks/{project}/{sample}.m6atm.benchmark.txt"
    params:
        extra=""
    conda:
        "../envs/m6atm.yaml"
    shell:
        "m6atm predict "
        "--bam {input.bam} "
        "--reference {input.reference} "
        "--config {input.config} "
        "--threads {threads} "
        "{params.extra} "
        "--output {output.predictions} "
        "2>{log} && touch {output.completion}"


rule m6atm_postprocess:
    input:
        predictions="{project}/results/m6atm/{sample}_predictions.completed",
        pred_file="{project}/results/m6atm/{sample}_predictions.tsv"
    output:
        directory="{project}/results/dataprep/{sample}_m6atm_dataprep")
    threads: 1
    resources:
        mem_mb = 1024 * 20
    priority: 10
    log:
        "logs/{project}/m6atm_postprocess/{sample}.log"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/m6atm_postprocess.py"
