rule differr:
    input:
        control_bam="results/alignments/{control}_filtered.bam",
        control_bai="results/alignments/{control}_filtered.bam.bai",
        native_bam="results/alignments/{native}_filtered.bam",
        native_bai="results/alignments/{native}_filtered.bam.bai",
        reference=config['reference']['transcriptome_fasta']
    output:
        "results/differr/{native}_{control}/{native}_{control}.differr.bed"
    params:
        prefix="{native}_{control}",
        extra=config['params']['differr']
    log:
        stdout="logs/differr/{native}_{control}.log"
    threads: config['threads']['differr']
    benchmark:
        "benchmarks/{native}_{control}.differr.benchmark.txt"
    container:
        "docker://btrspg/differr:latest"
    shell:
        "differr "
        " -a {input.control_bam} "
        " -b {input.native_bam} "
        " -r {input.reference} "
        " -o {output} "
        " {params.extra} "
        " -p {threads} 2>{log}"




rule differr_sampled:
    input:
        control_bam="results/alignments/{control}_filtered_{sample_size}_{n}.bam",
        control_bai="results/alignments/{control}_filtered_{sample_size}_{n}.bam.bai",
        native_bam="results/alignments/{native}_filtered_{sample_size}_{n}.bam",
        native_bai="results/alignments/{native}_filtered_{sample_size}_{n}.bam.bai",
        reference=config['reference']['transcriptome_fasta']
    output:
        "results/differr/{native}_{control}-{sample_size}-{n}/differr.bed"
    params:
        prefix="",
        extra=config['params']['differr']
    log:
        stdout="logs/differr/{native}_{control}-{sample_size}-{n}.log"
    threads: config['threads']['differr']
    benchmark:
        "benchmarks/{native}_{control}-{sample_size}-{n}.differr.benchmark.txt"
    container:
        "docker://btrspg/differr:latest"
    shell:
        "differr "
        " -a {input.control_bam} "
        " -b {input.native_bam} "
        " -r {input.reference} "
        " -o {output} "
        " {params.extra} "
        " -p {threads} 2>{log}"
