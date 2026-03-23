rule variants_bcftools:
    input:
        bam="{project}/results/alignments/{sample}.bam",
        reference=config["reference"]["transcriptome_fasta"],
    output:
        variants_vcf="{project}/results/variants/{sample}.bcf",
    container:
        get_container("bcftools")
    benchmark:
        "benchmarks/{project}/{sample}.variants_bcftools.benchmark.txt"
    log:
        "logs/{project}/qc/{sample}_variants_bcftools.log",
    priority: 20
    threads: get_threads("bcftools", 4)
    resources:
        mem_mb=1024 * 10,
        disk_mb=1024 * 10,
    shell:
        "bcftools mpileup --threads {threads} -f {input.reference} "
        "{input.bam} |bcftools call --threads {threads} -Ob -o {output.variants_vcf} -c 2>{log}"
