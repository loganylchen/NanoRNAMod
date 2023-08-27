rule link_fastq:
    input:
        fastq="data/{sample}/fastq/pass.fq.gz"
    output:
        fastq="results/fastq/{sample}.fq.gz",
    params:
        relative_path="../../data/{sample}/fastq/pass.fq.gz"
    log:
        "logs/link_fastq/{sample}.log"
    shell:
        "ln -s {params.relative_path} {output.fastq} && "
        "echo `date` > {log} "

rule read_assignment:
    input:
        read_assign_pkl="results/assembly/{sample}.lafiteread_assignment.pkl",
        fastq="results/fastq/{sample}.fq.gz"
    output:
        read_assignment_dir=temp(directory("results/read_assignment/{sample}_tmp")),
        mapping_list="results/read_assignment/{sample}.list",
    params:
        transcriptome_fasta=config['reference']['transcriptome_fasta'],
    log:
        "logs/read_assignment/{sample}.log",
        "logs/read_assignment/{sample}.err",
    benchmark:
        "benchmarks/{sample}.python_read_assignment.benchmark.txt"
    conda:
        "../envs/pyfastx.yaml"
    script:
        "../scripts/read_assignment.py"