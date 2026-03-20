import os
import yaml
import sys

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = open(snakemake.log[0], "w")

bam_file = snakemake.input.bam
reference = snakemake.input.reference
output_config = snakemake.output.config
modifications = snakemake.params.modifications
min_prob = snakemake.params.min_prob

config_dict = {
    'bam': bam_file,
    'reference': reference,
    'modifications': modifications,
    'min_prob': min_prob,
    'output_prefix': os.path.splitext(bam_file)[0]
}

with open(output_config, 'w') as f:
    yaml.dump(config_dict, f)

print(f"Configuration written to {output_config}")
print(f"Modifications: {modifications}")
print(f"Min probability: {min_prob}")
