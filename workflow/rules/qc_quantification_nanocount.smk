rule quant_nanocount:
    input:
        bam="{project}/results/alignments/{sample}.bam",
    output:
        count_tsv="{project}/results/quantification/{sample}.tx_counts.tsv",
    benchmark:
        "benchmarks/{project}/{sample}.quantification_nanocount.benchmark.txt"
    log:
        "logs/{project}/qc/{sample}_quantification_nanocount.log",
    container:
        "docker://btrspg/nanocount:latest"
    shell:
        "NanoCount -i {input.bam}  "
        "--extra_tx_info "
        "-o {output.count_tsv} > {log}"
