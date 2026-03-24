# py-baleen: CUDA-accelerated DTW + HMM modification detection
# https://github.com/loganylchen/py-baleen

rule pybaleen_run:
    input:
        native_bam="{project}/results/alignments/{native}.bam",
        native_fastq="{project}/results/fastq/{native}.fq.gz",
        native_blow5="{project}/results/blow5/{native}.blow5",
        control_bam="{project}/results/alignments/{control}.bam",
        control_fastq="{project}/results/fastq/{control}.fq.gz",
        control_blow5="{project}/results/blow5/{control}.blow5",
        ref=config["reference"]["transcriptome_fasta"],
    output:
        result=KEEP_OR_NOT(
            "{project}/results/pybaleen/{native}_{control}/site_results.tsv"
        ),
        pkl=KEEP_OR_NOT(
            "{project}/results/pybaleen/{native}_{control}/pipeline_results.pkl"
        ),
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
    priority: 5
    log:
        out="logs/{project}/pybaleen/N_{native}_C_{control}.log",
        err="logs/{project}/pybaleen/N_{native}_C_{control}.error",
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
        "{params.extra} "
        "1> {log.out} 2> {log.err}"
