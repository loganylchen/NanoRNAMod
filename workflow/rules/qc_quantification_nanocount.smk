rule quant_nanocount:
    input:
        bam="{project}/results/alignments/{sample}.bam",
    output:
        count_tsv="{project}/results/quantification/{sample}.tx_counts.tsv",
    benchmark:
        "benchmarks/{project}/{sample}.quantification_nanocount.benchmark.txt"
    log:
        stdout="logs/{project}/qc/{sample}_quantification_nanocount.log",
        stderr="logs/{project}/qc/{sample}_quantification_nanocount.err",
    container:
        "docker://btrspg/nanocount:latest"
    priority: 20
    resources:
        mem_mb=1024 * 10,
        disk_mb=1024 * 10,
    shell:
        "NanoCount -i {input.bam}  "
        "--extra_tx_info "
        "-o {output.count_tsv} 1> {log.stdout} 2>{log.stderr}"


rule qc_qualimap:
    input:
        bam="{project}/results/alignments/{sample}.splice.bam",
        bai="{project}/results/alignments/{sample}.splice.bam.bai",
    output:
        report="{project}/results/qc/{sample}/{sample}_rnaseq.pdf",
    params:
        output_dir=directory("{project}/results/qc/{sample}"),
        reference_gtf=config["reference"]["transcriptome_gtf"],
        output_file="{sample}_rnaseq.pdf",
    benchmark:
        "benchmarks/{project}/{sample}.qualimap_rnaseq.benchmark.txt"
    log:
        "logs/{project}/qc/{sample}_qualimap.log",
    conda:
        "../envs/qualimap.yaml"
    resources:
        mem="50G",
        mem_mb=1024 * 50,
        javaopt=" -Djava.awt.headless=true ",
    shell:
        "export JAVA_OPTS='{resources.javaopt}' && "
        "qualimap rnaseq -bam {input.bam} "
        "-gtf {params.reference_gtf} "
        "-outdir {params.output_dir} "
        "-outfile {params.output_file} "
        "-outformat PDF:HTML --java-mem-size={resources.mem}> {log}"
