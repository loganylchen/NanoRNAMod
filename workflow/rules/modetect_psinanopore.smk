rule psinanocompore_psi:
    input:
        control_bam="results/alignments/{control}.splice.bam",
        native_bam="results/alignments/{native}.splice.bam",
    output:
        output="results/psinanopore/{native}_{control}.psi_candidates.csv",
    log:
        "logs/psinanopore/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.psinanopore.benchmark.txt"
    container:
        "docker://btrspg/psinanopore:latest"
    params:
        reference=config['reference']['genome_fasta'],
        pvalue=config['params']['psinanopore'],
    threads: config['threads']['psinanopore']
    shell:
        "PsiDetect.R -f {input.native_bam} "
        "-g {input.control_bam} "
        "-k /opt/psinanopore/data/kmer_summary.csv "
        "-r {params.reference} "
        "-m {params.pvalue}  -o ~/Desktop/psi_candidates.csv"
