import sys
import os
import pysam
import gzip
import multiprocessing
import subprocess
import random

sys.stdout = open(snakemake.log[0], "w")

threads = snakemake.threads
sample_sizes = [int(i) for i in snakemake.params.sample_size]

sample_read_names = dict()

with pysam.AlignmentFile(snakemake.input.bam, "rb") as bam:
    for contig in bam.references:
        contig_len = bam.get_reference_length(contig)
        read_names = []
        for pc in bam.pileup(contig, int(contig_len/2)-2, int(contig_len/2)+2):
            read_names += pc.get_query_names()
        read_names = sorted(list(set(read_names)))
        for sample_size in sample_sizes:
            if len(read_names) > sample_size:
                sample_read_names[sample_size] = random.sample(read_names, sample_size)
            else:
                sample_read_names[sample_size] = read_names

param = list()
for read_name_sampled,output_bam in zip(snakemake.output.read_names,snakemake.output.bams):
    sample_size = int(os.path.basename(read_name_sampled).split('_')[-1].split('.')[0])
    with open(read_name_sampled, "w") as f:
        for read_name in set(sample_read_names[sample_size]):
            f.write(f"{read_name}\n")
    param.append([snakemake.input.bam, read_name_sampled, output_bam])

def samtools_view(input_bam, read_name_file, output_bam):
    cmd = f"samtools view -b -o {output_bam} -N {read_name_file} {input_bam} ; samtools index {output_bam}"
    p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr = p.communicate()
    print(f'''
    {cmd}
    stdout: {stdout.decode('utf8')}
    stderr: {stderr.decode('utf8')}
    ''')


with multiprocessing.Pool(threads) as pool:
    pool.starmap(samtools_view, param)

print('Done!')











