rule epinano_prep:
    input:
        sample_bam="results/alignments/{sample}_3.2.4_filtered.bam",
        sample_bai="results/alignments/{sample}_3.2.4_filtered.bam.bai",
        reference=config['reference']['transcriptome_fasta'],
        reference_dict=config['reference']['transcriptome_fasta'] + '.dict'
    output:
        per_site = temp("results/alignments/{sample}_3.2.4_filtered.plus_strand.per.site.csv"),
        sum_err=temp("results/alignments/{sample}_3.2.4_filtered.plus.sumErrOut.csv"),
        # kmer_5_site = "results/alignments/{sample}_filtered.plus_strand.per.site.5mer.csv",
        # dump_csv = "results/alignments/{sample}_filtered.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"
    params:
        extra=config['params']['epinano_dataprep'],
        prefix="results/alignments/{sample}_filtered"
    threads: config['threads']['epinano']
    log:
        "logs/epinano_prep/{sample}.log"
    benchmark:
        "benchmarks/{sample}.epinano_prep.benchmark.txt"
    container:
        "docker://btrspg/epinano:latest"
    script:
        "../scripts/epinano_prep.sh"






rule epinano:
    input:
        control="results/alignments/{control}_3.2.4_filtered.plus_strand.per.site.csv",
        native="results/alignments/{native}_3.2.4_filtered.plus_strand.per.site.csv",
    output:
        "results/epinano/{native}_{control}/epinano.delta-sum_err.prediction.csv"
    params:
        extra=config['params']['epinano'],
        prefix="results/epinano/{native}_{control}/epinano"
    threads: config['threads']['epinano']
    log:
        "logs/epinano/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.epinano.benchmark.txt"
    container:
        "docker://btrspg/epinano:latest"
    shell:
        "Rscript /opt/Epinano/Epinano_DiffErr.R "
        "-k {input.control} "
        "-w {input.native} "
        "-o {params.prefix} "
        "{params.extra} 2>{log}"



