import sys
import subprocess
from collections import defaultdict

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout


bam_sublist = defaultdict(list)

with open(snakemake.input.mapping_list, 'r') as mapping_list:
    i = 0
    for line in mapping_list:
        i += 1
        transcript, fastq, fasta, bam = line.strip().split('\t')
        if i % 1000 == 0:
            print(f'Working on {transcript}, {i}')
        bam_sublist[int(i/1000)].append(bam)
        # bam_list.write(f'{bam}\n')
        subprocess.Popen(
            f"minimap2 -t {snakemake.threads} {snakemake.params.extra} {fasta} {fastq} | samtools view -Sbh | samtools sort - -o {bam}",
            shell=True, stderr=sys.stdout, stdout=sys.stdout).communicate()
        subprocess.Popen(f'samtools index {bam}', shell=True, stderr=sys.stdout, stdout=sys.stdout).communicate()

with open(f"{snakemake.params.transcriptome_fasta}.fai", 'r') as fai, open(f"{snakemake.input.mapping_dir}/header.tmp",
                                                                           'w') as header:
    for line in fai:
        transcript, length, *_ = line.strip().split('\t')
        header.write(f'@SQ\tSN:{transcript}\tLN:{length}\n')

print('Merging')

with open(snakemake.output.bam_list, 'w') as bam_list, open(f'{snakemake.output.bam_list}.2round.tmp', 'w') as round2merge_list:
    for i,bam_sublist_key in enumerate(bam_sublist.keys()):
        print(f'As too many bams listing here, so we seperate them into different subbam, this is the {i+1} bam')
        with open(f'{snakemake.output.bam_list}.tmp', 'w') as tmp_bam_list:
            tmp_bam_list.write('\n'.join(bam_sublist[bam_sublist_key]))
            bam_list.write('\n'.join(bam_sublist[bam_sublist_key]))
        subprocess.Popen(
            f"samtools merge -b {snakemake.output.bam_list}.tmp -h {snakemake.input.mapping_dir}/header.tmp -o - | samtools sort - -o {snakemake.output.bam}.tmp{i}",
            shell=True, stderr=sys.stdout, stdout=sys.stdout).communicate()
        subprocess.Popen(f"samtools index {snakemake.output.bam}.tmp{i}", shell=True, stderr=sys.stdout,
                         stdout=sys.stdout).communicate()
        round2merge_list.write(f'{snakemake.output.bam}.tmp{i}\n')

if len(bam_sublist.keys()) > 1:
    subprocess.Popen(
        f"samtools merge -b {snakemake.output.bam_list}.2round.tmp -h {snakemake.input.mapping_dir}/header.tmp -o - | samtools sort - -o {snakemake.output.bam}",
        shell=True, stderr=sys.stdout, stdout=sys.stdout).communicate()
    subprocess.Popen(f"samtools index {snakemake.output.bam}", shell=True, stderr=sys.stdout,
                     stdout=sys.stdout).communicate()
    subprocess.Popen(f"rm {snakemake.output.bam_list}.2round.tmp {snakemake.input.mapping_dir}/header.tmp ",
                     shell=True, stderr=sys.stdout, stdout=sys.stdout).communicate()
else:
    subprocess.Popen(
        f"mv {snakemake.output.bam}.tmp{i} {snakemake.output.bam}",
        shell=True, stderr=sys.stdout, stdout=sys.stdout).communicate()
    subprocess.Popen(f"mv  {snakemake.output.bam}.tmp{i}.bai {snakemake.output.bam}.bai", shell=True, stderr=sys.stdout,
                     stdout=sys.stdout).communicate()
    subprocess.Popen(f"rm {snakemake.output.bam_list}.2round.tmp {snakemake.input.mapping_dir}/header.tmp ",
                     shell=True, stderr=sys.stdout, stdout=sys.stdout).communicate()


