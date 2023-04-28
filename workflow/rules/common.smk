import glob
import pandas as pd
import sys
from snakemake.utils import validate
from snakemake.logging import logger

validate(config, schema="../schemas/config.schema.yaml")

window_length1 = [3,4,5,6,7]
window_length2 = [8,9,10]
thresholds = [1+i/10 for i in range(4,20)]
peaks = [i/10 for i in range(1,20)]



# loading samples
samples = (
    pd.read_csv(config['samples'], sep="\t", dtype={"SampleName": str})
    .set_index("SampleName", drop=False)
    .sort_index()
)

validate(samples,schema='../schemas/samples.schema.yaml')

if len(samples['Condition'].unique()) != 2:
    logger.error('Samples should at least come from 2 conditions')
    sys.exit(1)
else:
    native_samples = list(samples.loc[samples['Condition']=='Native',:].index)
    control_samples = list(samples.loc[samples['Condition']=='Control',:].index)
    comparisons = [f'{ns}_{cs}' for ns in native_samples for cs in control_samples]


# def get_fastq():
#     dir_dict = samples.to_dict()['Directory']

def get_final_output():
    tools = [tool for tool in config['tools'] if config['tools'][tool]['activate']]
    final_output = expand("results/modifications/{comp}/{tool}.tsv.gz",comp=comparisons,tool=tools)
    final_output += expand("results/eventalign/{sample}_baleen_{w1}_{w2}_{threshold}_{peak}.tsv.bz2",sample=samples.index.values,w1=window_length1,w2=window_length2,threshold=thresholds,peak=peaks)
    return final_output



# def get_final_output():
#     final_output = expand(
#         "results/diffexp/{contrast}.diffexp.symbol.tsv",
#         contrast=config["diffexp"]["contrasts"],
#     )
#     final_output.append("results/deseq2/normcounts.symbol.tsv")
#     final_output.append("results/counts/all.symbol.tsv")
#     return final_output


    # expand("modifications/tombo_mod_lsc/{control}_{native}.lsc.tombo.stats",control=config['conditions'][
#             'control'],native=config['conditions']['native']),
#         expand("modifications/tombo_mod_msc/{control}_{native}.msc.tombo.stats",control=config['conditions'][
#             'control'],native=config['conditions']['native']),
#         expand("modifications/tombo_mod_msc/{control}_{native}.msc.tombo.per_read_stats",control=config['conditions'][
#             'control'],native=config['conditions']['native']),
#         expand("modifications/xpore/{control}_{native}/majority_direction_kmer_diffmod.table",control=
#         config['conditions']['control'],native=config['conditions']['native']),
#         expand("modifications/nanocompore/{control}_{native}",control=config['conditions']['control'],native=
#         config['conditions']['native']),
#         expand("modifications/m6anet/{native}",native=config['conditions']['native'])