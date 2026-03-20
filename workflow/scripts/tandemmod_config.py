import os
import yaml
import sys

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = open(snakemake.log[0], "w")

bam_file = snakemake.input.bam
reference = snakemake.input.reference
output_config = snakemake.output.config
model_type = snakemake.params.model_type
threshold = snakemake.params.threshold

config_dict = {
    'bam': bam_file,
    'reference': reference,
    'model_type': model_type,
    'threshold': threshold,
    'output_prefix': os.path.splitext(bam_file)[0]
}

with open(output_config, 'w') as f:
    yaml.dump(config_dict, f)

print(f"Configuration written to {output_config}")
print(f"Model type: {model_type}")
print(f"Threshold: {threshold}")
