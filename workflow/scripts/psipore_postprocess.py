import pandas as pd

input_file = snakemake.input[0]
output_file = snakemake.output[0]

df = pd.read_csv(input_file, sep='\t')
if 'transcript' in df.columns and 'position' in df.columns:
    df = df.sort_values(['transcript', 'position'])
df.to_csv(output_file, sep='\t', index=False)
