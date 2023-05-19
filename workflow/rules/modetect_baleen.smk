rule baleen_rmd:
    input:
        native_bam="results/alignments/{native}_filtered.bam",
        control_bam="results/alignments/{control}_filtered.bam",
        native_fast5_dir="data/{native}/fast5/workspace",
        control_fast5_dir="data/{control}/fast5/workspace",
        native_sequencing_summary="data/{native}/fast5/sequencing_summary.txt",
        control_sequencing_summary="data/{control}/fast5/sequencing_summary.txt",
        reference=config['reference']['transcriptome_fasta']
    output:
        outdir=directory("results/baleen/{native}_{control}"),
        completion="results/baleen/{native}_{control}/baleen_rmd.done"

    log:
        "logs/baleen/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.baleen_rmd.benchmark.txt"
    params:
        extra=config['params']['baleen']
    threads: config['threads']['baleen']
    container:
        "docker://btrspg/baleen:dev"
    shell:
        "baleenRMD "
        "--native-bam {input.native_bam} "
        "--control-bam {input.control_bam} "
        "--native-fast5-dir {input.native_fast5_dir} "
        "--control-fast5-dir {input.control_fast5_dir} " 
        "--native-fast5-filename {input.native_sequencing_summary} " 
        "--control-fast5-filename {input.control_sequencing_summary} "
        "--reference {input.reference} " 
        "--output-dir {output.outdir} "
        "--threads {threads} "
        "{params.extra} 2> {log} && echo `date` > {output.completion}"

rule baleen_dataprep:
    input:
        eventalign="results/eventalign/{sample}_nanocompore.tsv.gz",
        completion="results/eventalign/{sample}_nanocompore.tsv.completed"
    output:
        directory("results/dataprep/{sample}_baleen_dataprep"),
        "results/dataprep/{sample}_baleen_dataprep/dataprep/{sample}/data.index"
    params:
        label="{sample}"
    container:
        "docker://btrspg/baleen:dev"
    threads: config['threads']['baleen']
    log:
        "logs/baleen_dataprep/{sample}.log",
        "logs/baleen_dataprep/{sample}.error"
    script:
        "../scripts/baleen_dataprep.py"




