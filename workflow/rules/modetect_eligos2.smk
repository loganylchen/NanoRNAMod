rule eligos2:
    input:
        control_bam="results/alignments/{control}_filtered.bam",
        control_bai="results/alignments/{control}_filtered.bam.bai",
        native_bam="results/alignments/{native}_filtered.bam",
        native_bai="results/alignments/{native}_filtered.bam.bai",
        reference=config['reference']['transcriptome_fasta'],
        region="results/eligos2/{native}_{control}.bed"
    output:
        directory=directory("results/eligos2/{native}_{control}"),
        result="results/eligos2/{native}_{control}/{native}_filtered_vs_{control}_filtered_on_{native}_{control}_baseExt0.txt",
    params:
        prefix="{native}_{control}",
        extra=config['params']['eligos2']
    threads: config['threads']['eligos2']
    benchmark:
        "benchmarks/{native}_{control}.eligos2.benchmark.txt"
    container:
        "docker://piroonj/eligos2:latest"
    shell:
        "eligos2 pair_diff_mod "
        "-tbam {input.native_bam} "
        "-cbam {input.control_bam} "
        "-ref {input.reference} "
        "-t {threads} "
        "-reg {input.region} "
        "-o {output.directory} {params.extra}"



rule eligos2_sampled:
    input:
        control_bam="results/alignments/{control}_filtered_{sample_size}_{n}.bam",
        control_bai="results/alignments/{control}_filtered_{sample_size}_{n}.bam.bai",
        native_bam="results/alignments/{native}_filtered_{sample_size}_{n}.bam",
        native_bai="results/alignments/{native}_filtered_{sample_size}_{n}.bam.bai",
        reference=config['reference']['transcriptome_fasta'],
        region="results/eligos2/{native}_{control}-{sample_size}-{n}.bed"
    output:
        directory=directory("results/eligos2/{native}_{control}-{sample_size}-{n}"),
        result="results/eligos2/{native}_{control}-{sample_size}-{n}/{native}_filtered_{sample_size}_{n}_vs_{control}_filtered_{sample_size}_{n}_on_{native}_filtered_{sample_size}_{n}_{control}_filtered_{sample_size}_{n}_baseExt0.txt",
    params:
        prefix="{native}_{control}",
        extra=config['params']['eligos2']
    threads: config['threads']['eligos2']
    benchmark:
        "benchmarks/{native}_{control}-{sample_size}_{n}.eligos2.benchmark.txt"
    container:
        "docker://piroonj/eligos2:latest"
    shell:
        "eligos2 pair_diff_mod "
        "-tbam {input.native_bam} "
        "-cbam {input.control_bam} "
        "-ref {input.reference} "
        "-t {threads} "
        "-reg {input.region} "
        "-o {output} {params.extra}"


rule eligos2_prep_sampled:
    input:
        control_bam="results/alignments/{control}_filtered_{sample_size}_{n}.bam",
        control_bai="results/alignments/{control}_filtered_{sample_size}_{n}.bam.bai",
        native_bam="results/alignments/{native}_filtered_{sample_size}_{n}.bam",
        native_bai="results/alignments/{native}_filtered_{sample_size}_{n}.bam.bai",
    output:
        region=temp("results/eligos2/{native}_{control}-{sample_size}-{n}.bed")
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools bamtobed -bed12 -i {input.control_bam}  > {output.region}.tmp && "
        "bedtools bamtobed -bed12 -i {input.native_bam} >> {output.region}.tmp && "
        "bedtools sort -i {output.region}.tmp > {output.region}.tmp2 && "
        "bedtools merge -i {output.region}.tmp2 > {output.region} && "
        "rm {output.region}.tmp {output.region}.tmp2"

rule eligos2_prep:
    input:
        control_bam="results/alignments/{control}_filtered.bam",
        control_bai="results/alignments/{control}_filtered.bam.bai",
        native_bam="results/alignments/{native}_filtered.bam",
        native_bai="results/alignments/{native}_filtered.bam.bai",
    output:
        region=temp("results/eligos2/{native}_{control}.bed")
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools bamtobed -bed12 -i {input.control_bam}  > {output.region}.tmp && "
        "bedtools bamtobed -bed12 -i {input.native_bam} >> {output.region}.tmp && "
        "bedtools sort -i {output.region}.tmp > {output.region}.tmp2 && "
        "bedtools merge -i {output.region}.tmp2 > {output.region} && "
        "rm {output.region}.tmp {output.region}.tmp2"