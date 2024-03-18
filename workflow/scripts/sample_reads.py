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







# read_names=expand("results/alignments/{sample}_filtered_{sample_size}.txt", region=config['sample_size']),
# rule sample_reads:
#     input:
#         bam="results/alignments/{sample}_filtered.bam",
#         bai="results/alignments/{sample}_filtered.bam.bai",
#     output:
#         bams=expand("results/alignments/{sample}_filtered_{sample_size}.bam", region=config['sample_size']),
#         bais=expand("results/alignments/{sample}_filtered_{sample_size}.bam.bai", region=config['sample_size']),
#     log:
#         "logs/samplereads/{sample}.log"
#     conda:
#         "../envs/pysam.yaml"
#     threads: config['sample_reads']['threads']
#     params:
#         sample_size=config['sample_size']
#     script:
#         "../scripts/sample_reads.py"



transcriptome_fasta=snakemake.params.transcriptome_fasta
fastq_sequence=snakemake.input.fastq
read_assign_pkl=snakemake.input.read_assign_pkl
read_assignment_dir=snakemake.output.read_assignment_dir
mapping_list = snakemake.output.mapping_list

os.makedirs(read_assignment_dir,exist_ok=True)



with open(read_assign_pkl,'rb') as f:
    data = pickle.load(f)


fasta_sequences = pyfastx.Fasta(transcriptome_fasta)
fastq_sequence = pyfastx.Fastq(fastq_sequence)
with open(mapping_list,'w') as mapping_list_f:
    for transcript in fasta_sequences.keys():
        if transcript in data:
            mapping_list_f.write(f'{transcript}\t{read_assignment_dir}/{transcript}.fq.gz'
                                 f'\t{read_assignment_dir}/{transcript}.fasta'
                                 f'\t{read_assignment_dir}/{transcript}.bam\n')
            with gzip.open(f'{read_assignment_dir}/{transcript}.fq.gz', 'wt') as infastq, open(
                    f'{read_assignment_dir}/{transcript}.fasta', 'w') as infasta:
                infasta.write(f'>{transcript}\n')
                infasta.write(fasta_sequences[transcript].seq)
                for read in data[transcript]:
                    infastq.write(f'{fastq_sequence[read].raw}')




