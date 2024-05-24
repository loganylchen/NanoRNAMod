rule quant_nanocount:
    input:
        bam="results/alignments/{sample}.bam",
    output:
        count_tsv="results/quantification/{sample}.tx_counts.tsv"
    benchmark:
        "benchmarks/{sample}.quantification_nanocount.benchmark.txt"
    log:
        "logs/qc/{sample}_quantification_nanocount.log"
    container:
        "docker://btrspg/nanocount:latest"
    shell:
        "NanoCount -i {input.bam}  "
        "--extra_tx_info "
        "-o {output.count_tsv} > {log}"


