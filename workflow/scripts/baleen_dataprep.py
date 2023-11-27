
import sys
import os

from baleen.utils.squiggle.nanopolish import eventalign_index

sys.stderr = open(snakemake.log[1], "w")
sys.stdout = open(snakemake.log[0], "w")

eventalign_file = snakemake.input.eventalign
outdir=snakemake.output[0]


compress_types = ['bz2','gz']
file_type = os.path.basename(eventalign_file).split('.')[-1]
if file_type not in compress_types:
    file_type = None

os.makedirs(outdir,exist_ok=True)

print(eventalign_index(eventalign_file,outdir,chunksize=100000,compress=file_type,threads=snakemake.threads))







