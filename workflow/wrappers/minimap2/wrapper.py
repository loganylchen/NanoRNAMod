from snakemake.shell import shell

mapping_list=snakemake.input.mapping_list
mapping_dir=snakemake.input.mapping_dir

bam_list=snakemake.output.bam_list
output_bam=snakemake.output.bam
bai=snakemake.output.bai
threads=snakemake.threads



extra = snakemake.params.get("extra", "")
transcriptome_fasta=snakemake.params.transcriptome_fasta
log = snakemake.log_fmt_shell(stdout=False, stderr=True)



with open(mapping_list, 'r') as mapping_list_f, open(bam_list,'w') as bam_list_f:
    for line in mapping_list_f:
        transcript, fastq,fasta,bam = line.strip().split('\t')
        print(f'Working on {transcript}')
        bam_list_f.write(f'{bam}\n')
        shell(f"(minimap2 -t {threads} {extra} {fasta} {fastq} | samtools view -Sbh | samtools sort - -o {bam}) 2>> {log}")
        shell(f'samtools index {bam} 2>> {log}')

with open(f"{transcriptome_fasta}.fai", 'r') as fai, open(f"{mapping_dir}/header.tmp", 'w') as header:
    for line in fai:
        transcript, length, *_ = line.strip().split('\t')
        header.write(f'@SQ\tSN:{transcript}\tLN:{length}\n')

print('Merging')
shell(f"samtools merge -b {bam_list} -h {mapping_dir}/header.tmp | samtools sort - -o {output_bam} ")
shell(f"samtools index {output_bam}")
shell(f"rm {mapping_dir}/header.tmp")




