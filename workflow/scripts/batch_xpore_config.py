
import yaml
import sys

sys.stderr = open(snakemake.log[0], "w")



xpore_cf = {
    'notes': 'xpore analysis',
    'out': snakemake.params[0],
    'data': {
        'CASE': {f'rep{i+1}':j for i,j in enumerate(snakemake.input.native_dirs)},
        'CONTROL': {f'rep{i+1}':j for i,j in enumerate(snakemake.input.control_dirs)},
    }
}
with open(snakemake.output.conf,'w') as f:
    yaml.dump(xpore_cf,f)






