rule qualimap_rnaseq:
    input:
        bam="results/alignments/{sample}.splice.bam",
    output:
        output_dir = directory("results/qc/{sample}"),
        report = "results/qc/{sample}/{sample}_rnaseq.pdf"
    params:
        reference_gtf=config['reference']['transcriptome_gtf'],
        output_file="{sample}_rnaseq.pdf"
    benchmark:
        "benchmarks/{sample}.qualimap_rnaseq.benchmark.txt"
    log:
        "logs/qc/{sample}_qualimap.log"
    conda:
        "../envs/qualimap.yaml"
    shell:
        "qualimap rnaseq -b {input.bam} "
        "-gtf {params.reference_gtf} "
        "-outdir {output.output_dir} "
        "-outfile {params.output_file} > {log}"


