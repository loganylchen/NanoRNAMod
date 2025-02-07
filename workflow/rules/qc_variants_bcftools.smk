rule variants_bcftools:
    input:
        bam="{project}/results/alignments/{sample}.bam",
        reference=config["reference"]["transcriptome.fa"],
    output:
        variants_vcf="{project}/results/variants/{sample}.vcf",
    benchmark:
        "benchmarks/{project}/{sample}.variants_bcftools.benchmark.txt"
    log:
        "logs/{project}/qc/{sample}_variants_bcftools.log",
    conda:
        "../envs/minimap2.yaml"
    priority: 20
    resources:
        mem_mb=1024 * 10,
        disk_mb=1024 * 10,
    shell:
        "bcftools mpileup -f {input.reference} "
        "{input.bam} |bcftools call -Ov -o {output.variants_vcf} -c 2>{log}"
