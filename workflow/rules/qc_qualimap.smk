import snakemake.resources

rule qualimap_rnaseq:
    input:
        bam="results/alignments/{sample}.splice.bam",
    output:
        report = "results/qc/{sample}/{sample}_rnaseq.pdf"
    params:
        output_dir = directory("results/qc/{sample}"),
        reference_gtf=config['reference']['transcriptome_gtf'],
        output_file="{sample}_rnaseq.pdf"
    benchmark:
        "benchmarks/{sample}.qualimap_rnaseq.benchmark.txt"
    log:
        "logs/qc/{sample}_qualimap.log"
    conda:
        "../envs/qualimap.yaml"
    resources:
        mem='50G',
        javaopt="-Djava.io.tmpdir=results/tmp/"
    shell:
        "export JAVA_OPTS='{resources.javaopt}' && "
        "qualimap rnaseq -bam {input.bam} "
        "-gtf {params.reference_gtf} "
        "-outdir {params.output_dir} "
        "-outfile {params.output_file} "
        "-outformat PDF:HTML --java-mem-size={resources.mem}> {log}"


