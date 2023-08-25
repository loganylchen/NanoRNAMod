import sys
import os
import pyfastx
import pickle
import gzip


sys.stderr = open(snakemake.log[1], "w")
sys.stdout = open(snakemake.log[0], "w")



os.makedirs(snakemake.output[0],exist_ok=True)




fasta_sequences= pyfastx.Fasta('/home/logan/Codes/pybaleen/tmp/Homo_sapiens.GRCh38.cdna.all.fa')
fastq_sequence=pyfastx.Fastq('/home/logan/Codes/pybaleen/tmp/lafite/hek293tMettle3KO.pass.fq.gz')
with open('/home/logan/Codes/pybaleen/tmp/lafite/hek293tMettl3KORep1.lafiteread_assignment.pkl','rb') as f:
    data=pickle.load(f)
for transcript in track(fasta_sequences.keys(), total=len(fasta_sequences.keys())):
    if transcript in data:
        os.makedirs('/home/logan/Codes/pybaleen/tmp/lafite/4baleen', exist_ok=True)
        with gzip.open(f'/home/logan/Codes/pybaleen/tmp/lafite/4baleen/{transcript}.fq.gz', 'wt') as infastq, open(
                f'/home/logan/Codes/pybaleen/tmp/lafite/4baleen/{transcript}.fasta', 'w') as infasta:
            infasta.write(f'>{transcript}\n')
            infasta.write(fasta_sequences[transcript].seq)
            for read in data[transcript]:
                # infastq.write(f'{read}\n')
                infastq.write(f'{fastq_sequence[read].raw}')

control_paths = {
    'data':snakemake.input.control_dataprep.replace('data.index','data.nason'),
    'data_index':snakemake.input.control_dataprep
}

native_paths = {
    'data':snakemake.input.native_dataprep.replace('data.index','data.nason'),
    'data_index':snakemake.input.native_dataprep
}

gtf_file=snakemake.params.transcriptome_gtf

if gtf_file == '':
    gtf_file = None

params = {
    'padding': 2,
    'proba_method' : 'unify',
    'normalized' :True,
    'window_size' : 20,
    'coverage' : 0.9
}

result_dir = snakemake.output[0]
os.makedirs(result_dir,exist_ok=True)
molecule_json_file=molecule_modification_based_on_clustering(control_paths,native_paths,result_dir,params,threads=snakemake.threads,verbose=False)
transcriptome_2_genome(molecule_json_file,gtf_file,result_dir,threads=snakemake.threads,verbose=False)

with open(os.path.join(snakemake.output[0],'done.txt'),'w') as f:
    f.write('Done!!!\n')








