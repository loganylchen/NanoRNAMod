rule lafite_assembly:
    input:
        bam="results/alignments/{sample}.splice.bam",
        bai="results/alignments/{sample}.splice.bam.bai",
        polya="results/polya/{sample}.tsv.gz",
    output:
        assembly_gtf="results/assembly/{sample}.lafite.gtf",
        read_assignment="results/assembly/{sample}.lafiteread_assignment.pkl",
    params:
        annotation=config['reference']['transcriptome_gtf'],
        genome=config['reference']['genome_fasta'],
        extra=config['params']['lafite']
    threads: config['threads']['lafite']
    log:
        "logs/lafite_assembly/{sample}.log"
    benchmark:
        "benchmarks/{sample}.lafite_assembly.benchmark.txt"
    container:
        "docker://btrspg/lafite:dev"
    shell:
        "lafite -b {input.bam} "
        "-g {params.annotation} "
        "-f {params.genome} "
        "-o {output.assembly_gtf} "
        "-p {input.polya} "
        "-t {threads}  --read-assignment --assign-known {params.extra} 2>{log} "
