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


def get_final_output():
    # tools = [tool for tool in config['tools'] if config['tools'][tool]['activate']]
    # final_output = expand("results/modifications/{comp}/{tool}.tsv.gz",comp=comparisons,tool=tools)
    final_output = expand("results/assembly/{sample}.lafite.gtf",sample=list(samples.index))
    final_output += expand("results/polya/{sample}.tsv.gz",sample=list(samples.index))
    final_output += expand("results/alignments/{sample}.bam",sample=list(samples.index))
    final_output += expand("results/alignments/{sample}.realign.bam",sample=list(samples.index))
    final_output += expand("results/read_assignment/{sample}.list",sample=list(samples.index))
    return final_output



