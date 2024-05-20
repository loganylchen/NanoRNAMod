rule eligos2:
    input:
        control_bam="results/alignments/{control}_filtered.bam",
        control_bai="results/alignments/{control}_filtered.bam.bai",
        native_bam="results/alignments/{native}_filtered.bam",
        native_bai="results/alignments/{native}_filtered.bam.bai",
        reference=config['reference']['transcriptome_fasta']
    output:
        "results/eligos2/{native}_{control}"
    params:
        prefix="{native}_{control}",
        extra=config['params']['eligo2']
    threads: config['threads']['eligo2']
    log:
        stdout="logs/eligos2/{native}_{control}.log"

    benchmark:
        "benchmarks/{native}_{control}.eligos2.benchmark.txt"
    container:
        "docker://btrspg/eligos2:latest"
    shell:
        "eligos2 pair_diff_mod -tbam {input.native_bam} "
        "-cbam {input.control_bam} "
        "-ref {input.reference} "
        "-t {threads} "
        "-o {output} {params.extra}"



rule eligos2_sampled:
    input:
        control_bam="results/alignments/{control}_filtered_{sample_size}_{n}.bam",
        control_bai="results/alignments/{control}_filtered_{sample_size}_{n}.bam.bai",
        native_bam="results/alignments/{native}_filtered_{sample_size}_{n}.bam",
        native_bai="results/alignments/{native}_filtered_{sample_size}_{n}.bam.bai",
        reference=config['reference']['transcriptome_fasta']
    output:
        "results/eligos2/{native}_{control}-{sample_size}_{n}"
    params:
        prefix="{native}_{control}",
        extra=config['params']['eligo2']
    threads: config['threads']['eligo2']
    log:
        stdout="logs/eligos2/{native}_{control}-{sample_size}_{n}.log"
    benchmark:
        "benchmarks/{native}_{control}-{sample_size}_{n}.eligos2.benchmark.txt"
    container:
        "docker://btrspg/eligos2:latest"
    shell:
        "eligos2 pair_diff_mod -tbam {input.native_bam} "
        "-cbam {input.control_bam} "
        "-ref {input.reference} "
        "-t {threads} "
        "-o {output} {params.extra}"
