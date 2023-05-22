

import sys
import os
from baleen.workflows import molecule_modification, transcriptome_2_genome


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


params_dict = {
    'PCA_P4_HAAR_2': {
        'padding':4,
        'wavelet_level':2,
        'wavelet':'haar',
        'algo':'PCA',
        'proba_method':'unify',
        'normalized':True
    },
    'PCA_P4_HAAR_3': {
        'padding':4,
        'wavelet_level':3,
        'wavelet':'haar',
        'algo':'PCA',
        'proba_method':'unify',
        'normalized':True
    },
    'PCA_P4_DB2_3': {
        'padding':4,
        'wavelet_level':3,
        'wavelet':'db2',
        'algo':'PCA',
        'proba_method':'unify',
        'normalized':True
    },
    'PCA_P4_DB2_2': {
        'padding':4,
        'wavelet_level':2,
        'wavelet':'db2',
        'algo':'PCA',
        'proba_method':'unify',
        'normalized':True
    },
    'HBOS_P4_HAAR_2': {
        'padding':4,
        'wavelet_level':2,
        'wavelet':'haar',
        'algo':'HBOS',
        'proba_method':'unify',
        'normalized':True
    },
    'HBOS_P4_HAAR_3': {
        'padding':4,
        'wavelet_level':3,
        'wavelet':'haar',
        'algo':'HBOS',
        'proba_method':'unify',
        'normalized':True
    },
    'HBOS_P4_DB2_3': {
        'padding':4,
        'wavelet_level':3,
        'wavelet':'db2',
        'algo':'HBOS',
        'proba_method':'unify',
        'normalized':True
    },
    'HBOS_P4_DB2_2': {
        'padding':4,
        'wavelet_level':2,
        'wavelet':'db2',
        'algo':'HBOS',
        'proba_method':'unify',
        'normalized':True
    },
}


for key in params_dict:
    result_dir = os.path.join(snakemake.output[0],key)
    molecule_json_file=molecule_modification(control_paths,native_paths,result_dir,params_dict[key],threads=snakemake.threads,verbose=False)
    transcriptome_2_genome(molecule_json_file,gtf_file,result_dir,threads=snakemake.threads,verbose=False)

with open(os.path.join(snakemake.output[0],'done.txt'),'w') as f:
    f.write('Done!!!\n')








