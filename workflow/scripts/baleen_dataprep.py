
import yaml
import sys
from baleen.workflows import data_prep


sys.stderr = open(snakemake.log[0], "w")
sys.stdout = open(snakemake.log[0], "a")


print(snakemake.input.eventalign,snakemake.output[0], snakemake.params.label,snakemake.threads)
data_prep(snakemake.input.eventalign,snakemake.output[0], snakemake.params.label,threads=snakemake.threads,verbose=False)







