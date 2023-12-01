rule baleen_reads_sampling:
    input:
        bamfile="results/alignments/{sample}_filtered.bam",
        csifile="results/alignments/{sample}_filtered.bam.csi",
    output:
        bamfile="results/alignments/{sample}_filtered_sampling.bam",
        baifile="results/alignments/{sample}_filtered_sampling.bam.bai",
    log:
        "logs/baleen/{sample}_sampling.log"
    benchmark:
        "benchmarks/{sample}.baleen.sampling.txt",
    container:
        "docker://btrspg/baleen:dev",
    params:
        extra=config["params"]["baleen_bf"],
    shell:
        "baleen-BF "
        "--bamfile {input.bamfile} "
        "--output {output.bamfile} "
        "{params.extra} "

rule baleen_run:
    input:
        native_bam="results/alignments/{native}.realign.bam" if config['realign'] else "results/alignments/{native}_filtered.bam",
        native_bai="results/alignments/{native}.realign.bam.bai" if config['realign'] else "results/alignments/{native}_filtered.bam.bai",
        native_eventalign="results/eventalign/{native}_baleen.tsv.bz2",
        control_eventalign="results/eventalign/{control}_baleen.tsv.bz2",
    output:
        outdir=directory("results/baleen/{native}_{control}"),
        results="results/baleen/{native}_{control}/modification/result_genome.tsv"
    params:
        reference = config['reference']['transcriptome_fasta'],
        gtf = config['reference']['transcriptome_gtf'],
        extra = config['params']['baleen']
    container:
        "docker://btrspg/baleen:dev"
    threads: config['threads']['baleen']
    log:
        "logs/baleen/{native}_{control}.log"
    shell:
        "Baleen.py --control-eventalign {input.control_eventalign} "
        "--native-eventalign {input.native_eventalign} "
        "--native-bam {input.native_bam} "
        "--gtf {params.gtf} "
        "--ref-fasta {params.reference} {params.extra} "
        "--threads {threads} --output-dir {output.outdir} > {log}"

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
        eventalign="results/eventalign/{sample}_baleen.tsv.bz2",
        completion="results/eventalign/{sample}_baleen.tsv.completed",
        gtf="config['reference']['transcriptome_gtf']",
        reference="config['reference']['transcriptome_fasta']"
    output:
        outdir=directory("results/dataprep/{sample}_baleen_dataprep"),
        data="results/dataprep/{sample}_baleen_dataprep/data.nason"
    container:
        "docker://btrspg/baleen:dev"
    threads: config['threads']['baleen']
    log:
        out="logs/baleen_dataprep/{sample}.log",
        err="logs/baleen_dataprep/{sample}.error"
    shell:
        "Baleen.py dataprep "
        "--eventalign-file {input.eventalign} "
        "--gtf {input.gtf} "
        "--ref-fasta {input.reference} "
        "--output-dir {output.outdir} "
        "--threads {threads} 1>{log.out} 2>{log.err} "





rule baleen_test:
    input:
        native_dataprep="results/dataprep/{native}_baleen_dataprep/dataprep/{native}/data.index",
        control_dataprep="results/dataprep/{control}_baleen_dataprep/dataprep/{control}/data.index"
        
    output:
        directory("results/baleen/{native}_{control}"),
        "results/baleen/{native}_{control}/done.txt"
    params:
        transcriptome_gtf=config['reference']['transcriptome_gtf']
    container:
        "docker://btrspg/baleen:dev"
    threads: config['threads']['baleen']
    log:
        "logs/baleen/N_{native}_C_{control}.log",
        "logs/baleen/N_{native}_C_{control}.error"
    script:
        "../scripts/baleen_mod.py"


