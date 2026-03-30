rule epinano_prep:
    input:
        sample_bam="{project}/results/alignments/{sample}.bam",
        sample_bai="{project}/results/alignments/{sample}.bam.bai",
        reference=config["reference"]["transcriptome_fasta"],
        reference_fai=config["reference"]["transcriptome_fasta"] + ".fai",
    output:
        per_site=temp(
            "{project}/results/alignments/{sample}.plus_strand.per.site.var.csv"
        ),
        sum_err=temp(
            "{project}/results/alignments/{sample}.plus_strand.per.site.var.csv.per.site.var.sum_err.csv"
        ),
    params:
        extra=config["params"]["epinano_dataprep"],
    container:
        get_container("epinano")
    threads: get_threads("epinano", 4)
    resources:
        mem_mb = 1024 * 50
    priority: 10
    log:
        "logs/{project}/epinano_prep/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.epinano_prep.benchmark.txt"
    script:
        "../scripts/epinano_prep.sh"


rule epinano:
    input:
        control="{project}/results/alignments/{control}.plus_strand.per.site.var.csv",
        native="{project}/results/alignments/{native}.plus_strand.per.site.var.csv",
    output:
        results="{project}/results/epinano/{native}_{control}/epinano.delta-sum_err.prediction.csv",
        output_dir=KEEP_OR_NOT(
            directory("{project}/results/epinano/{native}_{control}")
        ),
    params:
        extra=config["params"]["epinano"],
        prefix="{project}/results/epinano/{native}_{control}/epinano",
    container:
        get_container("epinano")
    priority: 10
    threads: get_threads("epinano", 4)
    resources:
        mem_mb = 1024 * 50
    log:
        "logs/{project}/epinano/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.epinano.benchmark.txt"
    script:
        "../scripts/epinano_differr.sh"
