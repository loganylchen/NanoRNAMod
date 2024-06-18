import pandas as pd


df = pd.read_csv(snakemake.input[0],sep="\t",header=None)
df.columns = ["chrom","pos"] + snakemake.params.samples
df.to_csv(snakemake.output[0],sep="\t",index=False)