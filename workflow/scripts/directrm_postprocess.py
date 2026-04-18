import os
import sys
import pandas as pd

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout

predictions_file = snakemake.input.pred_file
output_dir = snakemake.output.directory

os.makedirs(output_dir, exist_ok=True)

output_file = f"{output_dir}/predictions.tsv"
if os.path.exists(predictions_file) and os.path.getsize(predictions_file) > 0:
    df = pd.read_csv(predictions_file, sep='\t')
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Processed DirectRM predictions: {len(df)} modifications detected")
    print(f"Modifications found: {df['modification'].unique().tolist() if 'modification' in df.columns else 'N/A'}")
    print(f"Output written to {output_file}")
else:
    print(f"Warning: Predictions file missing or empty: {predictions_file}")
    open(output_file, 'w').close()
