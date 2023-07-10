

import sys
import os
from baleen.workflows import molecule_modification_based_on_clustering, transcriptome_2_genome


sys.stderr = open(snakemake.log[1], "w")
sys.stdout = open(snakemake.log[0], "w")



os.makedirs(snakemake.output[0],exist_ok=True)

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








