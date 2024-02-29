import sys
import os

from baleen.fio.eventalign import Eventalign

sys.stderr = open(snakemake.log.err, "w")
sys.stdout = open(snakemake.log.out, "w")

eventalign_file = snakemake.input.eventalign
outdir = os.path.dirname(snakemake.output.data)
label = snakemake.params.label

threads = snakemake.threads

if not snakemake.params.use_mem:
    os.environ['JOBLIB_TEMP_FOLDER'] = f'{outdir}/tmp'


eventalign = Eventalign(eventalign_file, threads=threads)
eventalign.print_info()
eventalign.index(outdir,chunksize=100000)
eventalign.print_info()
