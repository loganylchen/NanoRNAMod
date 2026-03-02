import os
import sys
import pandas as pd

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout

predictions_file = snakemake.input.pred_file
output_dir = snakemake.output.directory

os.makedirs(output_dir, exist_ok=True)

if os.path.exists(predictions_file):
    df = pd.read_csv(predictions_file, sep='\t')
    output_file = f"{output_dir}/predictions.tsv"
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Processed RNANO predictions: {len(df)} modifications detected")
    if 'modification' in df.columns:
        print(f"Modifications found: {df['modification'].unique().tolist()}")
    print(f"Output written to {output_file}")
else:
    print(f"Warning: Predictions file not found: {predictions_file}")
    open(f"{output_dir}/predictions.tsv", 'w').close()
