import glob
import pandas as pd
import sys
from snakemake.utils import validate
from snakemake.logging import logger

validate(config, schema="../schemas/config.schema.yaml")


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


def get_nanocompore_list(sample_list):
    nanocompore_list = [f"results/nanocompore_eventalign_collapse/{sample}/{sample}_eventalign_collapse.tsv" for sample in sample_list]
    return ','.join(nanocompore_list)

iter_number = range(int(config['test_times']))

def get_final_output():
    tools = [tool for tool in config['tools'] if config['tools'][tool]['activate']]
    final_output = []
    final_output += expand("results/quantification/{sample}.tx_counts.tsv",sample=list(samples.index))
    # For some small dataset (on limited transcripts), sampling may be a good choice, while for other big datasets, it may not be necessary
    if config['sample']:
        if 'nanocompore' in tools:
            final_output += expand("results/nanocompore/{comp}-{sample_size}-{n}/nanocompore_results.tsv",comp=comparisons,
                sample_size=config['sample_size'],n=iter_number)
        if 'xpore' in tools:
            final_output += expand("results/xpore/{comp}-{sample_size}-{n}/xpore_results.tsv",comp=comparisons,sample_size=config[
            'sample_size'],n=iter_number)
        if 'baleen' in tools:
            final_output += expand('results/baleen/{comp}-{sample_size}-{n}/baleen_results.tsv',comp=comparisons,sample_size=config[
            'sample_size'],n=iter_number)
        if 'differr' in tools:
            final_output += expand("results/differr/{comp}-{sample_size}-{n}/differr_results.tsv",comp=comparisons,sample_size=config[
            'sample_size'],n=iter_number)
        if 'drummer' in tools:
            final_output += expand("results/drummer/{comp}-{sample_size}-{n}/",comp=comparisons,sample_size=config[
            'sample_size'],n=iter_number)
    else:
        if 'baleen' in tools:
            final_output += expand('results/baleen/{comp}/baleen_results.tsv',comp=comparisons)
        if 'nanocompore' in tools:
            final_output += expand("results/nanocompore/{comp}/nanocompore_results.tsv",comp=comparisons)
        if 'xpore' in tools:
            final_output += expand("results/xpore/{comp}/xpore_results.tsv",comp=comparisons)
        if 'differr' in tools:
            final_output += expand('results/differr/{comp}/differr_results.tsv',comp=comparisons)
        if 'drummer' in tools:
            final_output += expand("results/drummer/{comp}/",comp=comparisons)
    return final_output



