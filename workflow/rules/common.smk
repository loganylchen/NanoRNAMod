import glob
import pandas as pd
import sys
from snakemake.utils import validate
from snakemake.logging import logger

validate(config, schema="../schemas/config.schema.yaml")


# loading samples
samples = (
    pd.read_csv(config['samples'], sep="\t", dtype={"SampleName": str},comment='#')
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


def get_nanocompore_list(sample_list):
    nanocompore_list = [f"results/nanocompore_eventalign_collapse/{sample}/{sample}_eventalign_collapse.tsv" for sample in sample_list]
    return ','.join(nanocompore_list)



def get_final_output():
    tools = [tool for tool in config['tools'] if config['tools'][tool]['activate']]
    final_output = []
    final_output += expand("results/quantification/{sample}.tx_counts.tsv",sample=list(samples.index))
    # For some small dataset (on limited transcripts), sampling may be a good choice, while for other big datasets, it may not be necessary

    if 'baleen' in tools:
        final_output += expand('results/modifications/baleen/{comp}/baleen_results.tsv',comp=comparisons)
    if 'nanocompore' in tools:
        final_output += expand("results/modifications/nanocompore/{comp}/nanocompore_results.tsv",comp=comparisons)
    if 'xpore' in tools:
        final_output += expand("results/modifications/xpore/{comp}/xpore_results.tsv",comp=comparisons)
    if 'differr' in tools:
        final_output += expand('results/modifications/differr/{comp}/differr_results.tsv',comp=comparisons)
    if 'drummer' in tools:
        final_output += expand("results/modifications/drummer/{comp}/drummer_results.tsv",comp=comparisons)
    if 'eligos2' in tools:
        final_output += expand("results/modifications/eligos2/{comp}/eligos2_results.tsv",comp=comparisons)
    if 'epinano' in tools:
        final_output += expand("results/modifications/epinano/{comp}/epinano_results.tsv",comp=comparisons)
    return final_output



