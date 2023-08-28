import sys
import os
import pyfastx
import pickle
import gzip



sys.stdout = open(snakemake.log[0], "w")
sys.stderr = open(snakemake.log[1], "w")

transcriptome_fasta=snakemake.params.transcriptome_fasta
fastq_sequence=snakemake.input.fastq
read_assign_pkl=snakemake.input.read_assign_pkl
read_assignment_dir=snakemake.output.read_assignment_dir
mapping_list = snakemake.output.mapping_list

os.makedirs(read_assignment_dir,exist_ok=True)



with open(read_assign_pkl,'rb') as f:
    data = pickle.load(f)

with (pyfastx.Fasta(transcriptome_fasta) as fasta_sequences, pyfastx.Fastq(fastq_sequence) as fastq_sequence,
      open(mapping_list,'w') as mapping_list_f):
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



