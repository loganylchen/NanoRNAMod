import yaml
import sys

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = open(snakemake.log[0], "w")

bam_file = snakemake.input.bam
reference = snakemake.input.reference
output_config = snakemake.output.config
model = snakemake.params.model
gpu = snakemake.params.gpu

config_dict = {
    'bam': bam_file,
    'reference': reference,
    'model': model,
    'gpu': gpu,
    'output_prefix': bam_file.replace('.bam', '')
}

with open(output_config, 'w') as f:
    yaml.dump(config_dict, f)

print(f"Configuration written to {output_config}")
print(f"Model: {model}")
print(f"GPU: {gpu}")
