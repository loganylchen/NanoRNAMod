import sys
import os

from baleen.fio.eventalign import Eventalign

sys.stderr = open(snakemake.log.err, "w")
sys.stdout = open(snakemake.log.out, "w")

eventalign_file = snakemake.input.eventalign
outdir = os.path.dirname(snakemake.output.data)
label = snakemake.params.label

threads = snakemake.threads

eventalign = Eventalign(eventalign_file, threads=threads)
eventalign.print_info()
eventalign.index(outdir)
eventalign.print_info()
