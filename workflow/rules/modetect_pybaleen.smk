# py-baleen: CUDA-accelerated DTW + HMM modification detection
# https://github.com/loganylchen/py-baleen

rule pybaleen_run:
    input:
        native_bam="{project}/results/alignments/{native}.bam",
        native_bai="{project}/results/alignments/{native}.bam.bai",
        native_fastq="{project}/results/fastq/{native}.fq.gz",
        native_blow5="{project}/results/blow5/{native}.blow5",
        native_blow5_index="{project}/results/blow5/{native}.blow5.idx",
        native_fastq_index=multiext(
            "{project}/results/fastq/{native}.fq.gz",
            ".index",
            ".index.fai",
            ".index.gzi",
            ".index.readdb",
        ),
        control_bam="{project}/results/alignments/{control}.bam",
        control_bai="{project}/results/alignments/{control}.bam.bai",
        control_fastq="{project}/results/fastq/{control}.fq.gz",
        control_blow5="{project}/results/blow5/{control}.blow5",
        control_blow5_index="{project}/results/blow5/{control}.blow5.idx",
        control_fastq_index=multiext(
            "{project}/results/fastq/{control}.fq.gz",
            ".index",
            ".index.fai",
            ".index.gzi",
            ".index.readdb",
        ),
        ref=config["reference"]["transcriptome_fasta"],
    output:
        result=temp("{project}/results/pybaleen/{native}_{control}/site_results.tsv"),
        read_bam=KEEP_OR_NOT("{project}/results/pybaleen/{native}_{control}/read_results.bam"),
        read_bai=KEEP_OR_NOT("{project}/results/pybaleen/{native}_{control}/read_results.bam.bai"),
    params:
        output_dir="{project}/results/pybaleen/{native}_{control}/",
        extra=config["params"].get("pybaleen", ""),
    container:
        get_container("pybaleen")
    benchmark:
        "benchmarks/{project}/{native}_{control}.pybaleen.txt"
    threads: get_threads("pybaleen", 4)
    resources:
        mem_mb=1024 * 650,
        gpu=1,
    priority: 50
    log:
        "logs/{project}/pybaleen/N_{native}_C_{control}.log",
    shell:
        "baleen run "
        "--native-bam {input.native_bam} "
        "--native-fastq {input.native_fastq} "
        "--native-blow5 {input.native_blow5} "
        "--ivt-bam {input.control_bam} "
        "--ivt-fastq {input.control_fastq} "
        "--ivt-blow5 {input.control_blow5} "
        "--ref {input.ref} "
        "-o {params.output_dir} "
        "--threads {threads} "
        "{params.extra} "
        "> {log} 2>&1"
