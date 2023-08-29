import sys
import subprocess
sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout





with open(snakemake.input.mapping_list, 'r') as mapping_list, open(snakemake.output.bam_list,'w') as bam_list:
   for line in mapping_list:
       transcript, fastq,fasta,bam = line.strip().split('\t')
       print(f'Working on {transcript}')
       bam_list.write(f'{bam}\n')
       subprocess.Popen(f"minimap2 -t {snakemake.threads} {snakemake.params.extra} {fasta} {fastq} | samtools view -Sbh | samtools sort - -o {bam}",shell=True,stderr=sys.stdout,stdout=sys.stdout).communicate()
       subprocess.Popen(f'samtools index {bam}',shell=True,stderr=sys.stdout,stdout=sys.stdout).communicate()

with open(f"{snakemake.params.transcriptome_fasta}.fai", 'r') as fai, open(f"{snakemake.input.mapping_dir}/header.tmp", 'w') as header:
   for line in fai:
       transcript, length, *_ = line.strip().split('\t')
       header.write(f'@SQ\tSN:{transcript}\tLN:{length}\n')

print('Merging')
subprocess.Popen(f"samtools merge -b {snakemake.output.bam_list} -h {snakemake.input.mapping_dir}/header.tmp -o {snakemake.input.mapping_dir}/header.tmp.bam",shell=True,stderr=sys.stdout,stdout=sys.stdout).communicate()
subprocess.Popen(f" samtools sort {snakemake.input.mapping_dir}/header.tmp.bam -o {snakemake.output.bam} ")
subprocess.Popen(f"samtools index {snakemake.output.bam}",shell=True,stderr=sys.stdout,stdout=sys.stdout).communicate()
subprocess.Popen(f"rm {snakemake.input.mapping_dir}/header.tmp {snakemake.input.mapping_dir}/header.tmp.bam",shell=True,stderr=sys.stdout,stdout=sys.stdout).communicate()





