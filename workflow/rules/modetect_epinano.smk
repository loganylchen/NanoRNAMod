rule epinano_prep:
    input:
        sample_bam="{project}/results/alignments/{sample}_filtered.bam",
        sample_bai="{project}/results/alignments/{sample}_filtered.bam.bai",
        reference=config["reference"]["transcriptome_fasta"],
        reference_dict=config["reference"]["transcriptome_fasta"] + ".dict",
    output:
        per_site=temp(
            "{project}/results/alignments/{sample}_filtered.plus_strand.per.site.csv"
        ),
        # sum_err=temp(
        #     "{project}/results/alignments/{sample}_filtered.plus.sumErrOut.csv"
        # ),
        # kmer_5_site = "results/alignments/{sample}_filtered.plus_strand.per.site.5mer.csv",
        # dump_csv = "results/alignments/{sample}_filtered.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"
    params:
        extra=config["params"]["epinano_dataprep"],
        prefix="{project}/results/alignments/{sample}_filtered",
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
    container:
        get_container("epinano")
    script:
        "../scripts/epinano_prep.sh"


rule epinano:
    input:
        control="{project}/results/alignments/{control}_filtered.plus_strand.per.site.csv",
        native="{project}/results/alignments/{native}_filtered.plus_strand.per.site.csv",
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
        stdout="logs/{project}/epinano/{native}_{control}.log",
        stderr="logs/{project}/epinano/{native}_{control}.err",
    benchmark:
        "benchmarks/{project}/{native}_{control}.epinano.benchmark.txt"
    container:
        get_container("epinano")
    shell:
        "Rscript /opt/epinano/Epinano_DiffErr.R "
        "-t {threads} "
        "-k {input.control} "
        "-w {input.native} "
        "-o {params.prefix} "
        "{params.extra} 1>{log.stdout} 2>{log.stderr}"
