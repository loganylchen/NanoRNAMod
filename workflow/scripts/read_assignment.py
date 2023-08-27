import sys
import os
import pyfastx
import pickle
import gzip
import json


sys.stdout = open(snakemake.log[0], "w")
sys.stderr = open(snakemake.log[1], "w")

transcriptome_fasta=snakemake.params.transcriptome_fasta
fastq_sequence=snakemake.input.fastq
read_assign_pkl=snakemake.input.read_assign_pkl
read_assignment_dir=snakemake.output.read_assignment_dir
mapping_json = snakemake.output.mapping_json

os.makedirs(read_assignment_dir,exist_ok=True)

mapping_dict=dict()

with open(read_assign_pkl,'rb') as f:
    data = pickle.load(f)

with pyfastx.Fasta(transcriptome_fasta) as fasta_sequences, pyfastx.Fastq(fastq_sequence) as fastq_sequence:
    for transcript in fasta_sequences.keys():
        if transcript in data:
            mapping_dict[transcript] = {
                'fastq':f'{read_assignment_dir}/{transcript}.fq.gz',
                'reference': f'{read_assignment_dir}/{transcript}.fasta'
            }
            with gzip.open(f'{read_assignment_dir}/{transcript}.fq.gz', 'wt') as infastq, open(
                    f'{read_assignment_dir}/{transcript}.fasta', 'w') as infasta:
                infasta.write(f'>{transcript}\n')
                infasta.write(fasta_sequences[transcript].seq)
                for read in data[transcript]:
                    infastq.write(f'{fastq_sequence[read].raw}')

with open(mapping_json,'w') as f:
    json.dump(mapping_dict,f)
