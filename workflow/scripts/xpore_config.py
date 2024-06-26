
import yaml
import sys

sys.stderr = open(snakemake.log[0], "w")

xpore_cf = {
    'notes': 'xpore analysis',
    'out': snakemake.params[0],
    'data': {
        'CASE': {'rep1': snakemake.input.native_dir},
        'CONTROL': {'rep1': snakemake.input.control_dir}
    }

}
with open(snakemake.output.conf,'w') as f:
    yaml.dump(xpore_cf,f)





