

import sys
import os


from baleen.algo.mod import modcall_molecule_joblib,modcall_transcript_joblib
from baleen.fio.eventalign import Eventalign
from baleen.utils.bed import readBED



sys.stderr = open(snakemake.log.err, "w")
sys.stdout = open(snakemake.log.out, "w")

threads=snakemake.threads


params = {
    'padding': 1,
    'proba_method': 'unify',
    'coverage': 0.8,
    'min_depth': 15,
    'determine_proba': 0.9,
    'dtw_normalization': None,
    'dedi_method': 'umap',
    'dedi_component_n': 2,
    'gmm_component_n': 2,
    'decide_method': 'gmm',
}


sample = 100
bedfile = snakemake.params.bedfile
outdir = snakemake.output.outdir
re_index = False


if not snakemake.params.use_mem:
    os.environ['JOBLIB_TEMP_FOLDER'] = f'{outdir}/tmp'

if (bedfile is not None) and (os.path.exists(bedfile)):
    target_regions = readBED(bedfile)
else:
    target_regions = None

native_eventalign = Eventalign(snakemake.input.native_eventalign, threads=threads)
native_eventalign.print_info()
native_eventalign.index(os.path.dirname(snakemake.input.native_eventalign_index), reindex=re_index)
native_eventalign.print_info()



control_eventalign = Eventalign(snakemake.input.control_eventalign, threads=threads)
control_eventalign.print_info()
control_eventalign.index(os.path.dirname(snakemake.input.control_eventalign_index), reindex=re_index)
control_eventalign.print_info()


modcall_molecule_joblib(native_eventalign, control_eventalign, target_regions, threads, sample, params, f'{outdir}/modcall_sm')
modcall_transcript_joblib(f'{outdir}/modcall_sm', outdir, threads)














