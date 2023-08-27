rule minimap2_transcriptome_align:
    input:
        fastq="results/fastq/{sample}.fq.gz",
    output:
        bam="results/alignments/{sample}.bam",
        csi="results/alignments/{sample}.bam.csi",
        bai="results/alignments/{sample}.bam.bai",
    log:
        "logs/minimap2_transcriptome_alignment/{sample}.log"
    benchmark:
        "benchmarks/{sample}.minimap2_transcriptome_alignment.benchmark.txt"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config["params"]["minimap2_transcriptome"],
        reference=config['reference']['transcriptome_fasta'],
    threads: config['threads']['minimap2']
    shell:
        "minimap2 -t {threads} {params.extra} {params.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} && "
        "samtools index {output.bam}"

rule minimap2_list_align:
    input:
        mapping_list="results/read_assignment/{sample}.list",
        mapping_dir="results/read_assignment/{sample}_tmp",
        sam_header="results/read_assignment/{sample}.samheader.txt"
    output:
        bam_list="results/read_assignment/{sample}.bamlist",
        bam="results/alignments/{sample}.realign.bam",
        bai="results/alignments/{sample}.realign.bam.bai",
    log:
        "logs/minimap2_list_alignment/{sample}.log"
    benchmark:
        "benchmarks/{sample}.minimap2_list_alignment.benchmark.txt"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config["params"]["minimap2_transcriptome"],
        transcriptome_fasta=config['reference']['transcriptome_fasta']
    threads: config['threads']['minimap2']
    shell:
        """
for line in $(cat {input.mapping_list})
do
    trnascript=$(echo $line|awk -F '\t' 'print $1');
    fastq=$(echo $line|awk -F '\t' 'print $2');
    fasta=$(echo $line|awk -F '\t' 'print $3');
    bam=$(echo $line|awk -F '\t' 'print $4');
    echo "Working on $transcript" >> {log} ;
    minimap2 -t {threads} {params.extra} $fasta $fastq | samtools view -Sbh | samtools sort - -o $bam;
    samtools index $bam
done
awk -F "\t" "print $4" {input.mapping_list} > {output.bam_list}
awk -F "\t" 'print @SQ\tSN:$1\tLN:$2' {params.transcriptome_fasta}.fai > {input.mapping_dir}/header.tmp
samtools merge -b {output.bam_list} -h {input.mapping_dir}/header.tmp | samtools sort - -o {output.bam} ;
samtools index ${output.bam}
rm {input.mapping_dir}/header.tmp
        """

rule minimap2_genome_align:
    input:
        fastq="results/fastq/{sample}.fq.gz",
    output:
        bam="results/alignments/{sample}.splice.bam",
        csi="results/alignments/{sample}.splice.bam.csi",
        bai="results/alignments/{sample}.splice.bam.bai"
    log:
        "logs/minimap2_genome_alignment/{sample}.log"
    benchmark:
        "benchmarks/{sample}.minimap2_genome_alignment.benchmark.txt"
    conda:
        "../envs/minimap2.yaml"
    params:
        extra=config["params"]["minimap2_genome"],
        reference=config['reference']['genome_fasta'],
    threads: config['threads']['minimap2']
    shell:
        "minimap2 -t {threads} {params.extra} {params.reference} {input.fastq} 2>> {log} "
        "| samtools view -Sbh "
        "| samtools sort - -o {output.bam} --write-index 2>>{log} && "
        "samtools index {output.bam}"
