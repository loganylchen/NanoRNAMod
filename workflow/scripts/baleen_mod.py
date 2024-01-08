

import sys
import os

from baleen.fio.eventalign import Eventalign
from baleen.utils.bed import readBED
from baleen.algo.mod import modcall_transcript_joblib

sys.stderr = open(snakemake.log.err, "w")
sys.stdout = open(snakemake.log.out, "w")

threads=snakemake.threads




native_eventalign_file= snakemake.input.native_eventalign
native_eventalign_index_file = snakemake.input.native_eventalign_index
control_eventalign_file= snakemake.input.control_eventalign
control_eventalign_index_file = snakemake.input.control_eventalign_index

native_eventalign_index_dir=os.path.dirname(native_eventalign_index_file)
control_eventalign_index_dir=os.path.dirname(control_eventalign_index_file)

bedfile=snakemake.params.bedfile
outdir = os.path.dirname(snakemake.output.result)


sample=100
params=dict(
    padding=1,
    proba_method='unify',
    coverage=0.95,
    min_depth=30,
    dtw_normalization='luscinia',
    dedi_method='tsne',
    dedi_component_n=2,
    gmm_component_n=1,
    cut_end=5
)


if os.path.exists(bedfile):
    target_regions = readBED(bedfile)
else:
    target_regions = dict()

native_eventalign = Eventalign(native_eventalign_file,threads=threads,force_no_check=True)
control_eventalign = Eventalign(control_eventalign_file,threads=threads,force_no_check=True)

native_eventalign.index(native_eventalign_index_dir)
control_eventalign.index(control_eventalign_index_dir)


modcall_transcript_joblib(native_eventalign,control_eventalign, target_regions,threads,sample,params,outdir)


with open(snakemake.output.result,'w') as f:
    f.write('Done!!!\n')








