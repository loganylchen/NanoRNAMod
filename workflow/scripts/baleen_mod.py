

import sys
import os
import json

from baleen.algo.mod import modcall_molecule_joblib,modcall_transcript_joblib
from baleen.fio.eventalign import EventalignIndex
from baleen.utils.bed import readBED



sys.stderr = open(snakemake.log.err, "w")
sys.stdout = open(snakemake.log.out, "w")

threads=snakemake.threads


params = {
    'padding': int(snakemake.params.padding),
    'proba_method': 'unify',
    'coverage': float(snakemake.params.coverage),
    'min_depth': 15,
    'proba_threshold': 0.9,
    'dtw_normalization': None if snakemake.params.dtw_normalization == 'None' else snakemake.params.dtw_normalization,
    'dedi_method': 'umap',
    'dedi_component_n': 2,
    'gmm_component_n': int(snakemake.params.gmm_component_n),
    'decide_method': 'gmm',
    'cut_end': 5,
    'segment_depth': 15,
    'sample': int(snakemake.params.sample),
}

sample = params['sample']

bedfile = snakemake.params.bedfile
outdir = os.path.dirname(snakemake.output.result)

with open(f'{outdir}/modcall_params.json', 'w') as f:
    json.dump(params, f)

if not snakemake.params.use_mem:
    os.environ['JOBLIB_TEMP_FOLDER'] = f'{outdir}/tmp'

if (bedfile is not None) and (os.path.exists(bedfile)):
    target_regions = readBED(bedfile)
else:
    target_regions = None


native_eventalign_index = EventalignIndex(os.path.dirname(snakemake.input.native_eventalign_index))
control_eventalign_index = EventalignIndex(os.path.dirname(snakemake.input.control_eventalign_index))



modcall_molecule_joblib(native_eventalign_index,control_eventalign_index, target_regions,threads,sample,params,f'{outdir}/modcall_sm')
modcall_transcript_joblib(f'{outdir}/modcall_sm', outdir, threads,params['proba_threshold'])














