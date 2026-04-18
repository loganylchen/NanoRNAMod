import os
import sys
import pandas as pd
import shutil

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout

predictions_file = snakemake.input.pred_file
output_dir = snakemake.output.directory

os.makedirs(output_dir, exist_ok=True)

output_file = f"{output_dir}/predictions.tsv"
if os.path.exists(predictions_file) and os.path.getsize(predictions_file) > 0:
    df = pd.read_csv(predictions_file, sep='\t')
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Processed TandemMod predictions: {len(df)} modifications detected")
    print(f"Output written to {output_file}")
else:
    # Either missing (shouldn't happen — rule touches sentinel on success)
    # or empty (tool succeeded but produced no results). Emit empty sentinel.
    print(f"Warning: Predictions file missing or empty: {predictions_file}")
    open(output_file, 'w').close()
