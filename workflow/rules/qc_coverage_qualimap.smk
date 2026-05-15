
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
    container:
        get_container("qualimap")
    benchmark:
        "benchmarks/{project}/{sample}.qualimap_rnaseq.benchmark.txt"
    log:
        "logs/{project}/qc/{sample}_qualimap.log",
    priority: 50
    resources:
        javaopt=" -Djava.awt.headless=true ",
        mem=50,  # Java heap size in GB for qualimap
    shell:
        "export JAVA_OPTS='{resources.javaopt}' && "
        "qualimap rnaseq -bam {input.bam} "
        "-gtf {params.reference_gtf} "
        "-outdir {params.output_dir} "
        "-outfile {params.output_file} "
        "-outformat PDF:HTML --java-mem-size={resources.mem}G > {log}"
