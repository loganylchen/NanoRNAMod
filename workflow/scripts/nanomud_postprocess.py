import pandas as pd

input_file = snakemake.input.predictions
output_file = snakemake.output[0]

df = pd.read_csv(input_file, sep='\t')
# NanoMUD reports a 'modification' column (e.g. 'Psi', 'm1Psi')
if 'transcript' in df.columns and 'position' in df.columns:
    sort_cols = ['transcript', 'position']
    if 'modification' in df.columns:
        sort_cols.append('modification')
    df = df.sort_values(sort_cols)
df.to_csv(output_file, sep='\t', index=False)
