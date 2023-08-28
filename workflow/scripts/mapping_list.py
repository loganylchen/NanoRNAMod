import sys
import os
import subprocess



sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout





with open(snakemake.input.mapping_list, 'r') as mapping_list, open(snakemake.output.bam_list,'w') as bam_list:
    for line in mapping_list:
        transcript, fastq,fasta,bam = line.strip().split('\t')
        print(f'Working on {transcript}')
        bam_list.write(f'{bam}\n')
        subprocess.check_call(f"minimap2 -t {snakemake.threads} {snakemake.params.extra} {fasta} {fastq} | samtools view -Sbh | samtools sort - -o {bam}",shell=True)
        subprocess.check_call(f'samtools index {bam}',shell=True)

with open(f"{snakemake.params.transcriptome_fasta}.fai", 'r') as fai, open(f"{snakemake.input.mapping_dir}/header.tmp", 'w') as header:
    for line in fai:
        transcript, length, *_ = line.strip().split('\t')
        header.write(f'@SQ\tSN:{transcript}\tLN:{length}\n')

print('Merging')
subprocess.check_call(f"samtools merge -b {snakemake.output.bam_list} -h {snakemake.input.mapping_dir}/header.tmp | samtools sort - -o {snakemake.output.bam} ",shell=True)
subprocess.check_call(f"samtools index {snakemake.output.bam}",shell=True)
subprocess.check_call(f"rm {snakemake.input.mapping_dir}/header.tmp",shell=True)




