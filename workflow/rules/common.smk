import glob
import pandas as pd
import sys
from snakemake.utils import validate
from snakemake.logging import logger

validate(config, schema="../schemas/config.schema.yaml")

wildcard_constraints:
    sample="[\da-zA-Z\.]+",
    control="[\da-zA-Z\.]+",
    native="[\da-zA-Z\.]+"



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
    control_list= '-'.join(control_samples)
    native_list = '-'.join(native_samples)

def get_nanocompore_list(sample_list):
    nanocompore_list = [f"results/nanocompore_eventalign_collapse/{sample}/{sample}_eventalign_collapse.tsv" for sample in sample_list]
    return ','.join(nanocompore_list)



def get_final_output():
    tools = [tool for tool in config['tools'] if config['tools'][tool]['activate']]
    final_output = []
    final_output += expand("results/dataprep/{sample}_baleen_dataprep/eventalign.index",sample=list(samples.index))
    # final_output += expand("results/modifications/{comp}/{tool}.tsv.gz",comp=comparisons,tool=tools)
    # final_output += expand("results/assembly/{sample}.lafite.gtf",sample=list(samples.index))
    # final_output += expand("results/polya/{sample}.tsv.gz",sample=list(samples.index))
    # final_output += expand("results/alignments/{sample}.bam",sample=list(samples.index))
    # final_output += expand("results/alignments/{sample}.realign.bam",sample=list(samples.index))
    # final_output += expand("results/read_assignment/{sample}.list",sample=list(samples.index))
    # final_output += expand("results/eventalign/{sample}_baleen.tsv.bz2",sample=list(samples.index))
    # final_output += expand("results/qc/{sample}/{sample}_rnaseq.pdf",sample=list(samples.index))
    # final_output += expand("results/quantification/{sample}.tx_counts.tsv",sample=list(samples.index))
    # final_output += [f"results/xpore/Groups/{native_list}_{control_list}.xpore_config.yaml"]
    # final_output += [f"results/xpore/Groups/{native_list}_{control_list}/majority_direction_kmer_diffmod.table",f"results/nanocompore/Group_{native_list}_{control_list}"]
    return final_output



