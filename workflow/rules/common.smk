import glob
import pandas as pd
import sys
from snakemake.utils import validate
from snakemake.logging import logger

validate(config, schema="../schemas/config.schema.yaml")

wildcard_constraints:
    sample="[\da-zA-Z\.]+",
    control="[\da-zA-Z\.]+",
    native="[\da-zA-Z\.]+",
    sample_size="[\d]+"



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

    # final_output += expand('results/baleen/{native}_{control}/done.txt',native=native_samples,control=control_samples)
    # final_output += expand("results/modifications/{comp}/{tool}.tsv.gz",comp=comparisons,tool=tools)
    # final_output += expand("results/assembly/{sample}.lafite.gtf",sample=list(samples.index))
    # final_output += expand("results/polya/{sample}.tsv.gz",sample=list(samples.index))
    # final_output += expand("results/alignments/{sample}.bam",sample=list(samples.index))
    # final_output += expand("results/alignments/{sample}.realign.bam",sample=list(samples.index))
    # final_output += expand("results/read_assignment/{sample}.list",sample=list(samples.index))
    # final_output += expand("results/eventalign/{sample}_baleen.tsv.bz2",sample=list(samples.index))
    # final_output += expand("results/qc/{sample}/{sample}_rnaseq.pdf",sample=list(samples.index))
    # final_output += expand("results/quantification/{sample}.tx_counts.tsv",sample=list(samples.index))
    #
    if config['sample']:
        final_output += expand("results/alignments/{sample}_filtered_{sample_size}.bam",sample=list(samples.index) ,sample_size=config[
            'sample_size']),
        if 'nanocompore' in tools:
            final_output += expand("results/nanocompore/{comp}-{sample_size}/{comp}_{sample_size}nanocompore_results.tsv",comp=comparisons,
                sample_size=config[
            'sample_size'])
        if 'xpore' in tools:
            final_output += expand("results/xpore/{comp}-{sample_size}/majority_direction_kmer_diffmod.table",comp=comparisons,sample_size=config[
            'sample_size'])
        if 'm6anet' in tools:
            final_output += expand("results/m6anet/{sample}-{sample_size}/data.site_proba.csv",sample=list(samples.index),sample_size=config[
            'sample_size'])
        if 'psinanopore' in tools:
            final_output += expand("results/psinanopore/{comp}-{sample_size}.psi_candidates.csv",comp=comparisons,sample_size=config[
            'sample_size'])
        if 'baleen' in tools:
            final_output += expand('results/baleen/{comp}-{sample_size}/transcripts.csv',comp=comparisons,sample_size=config[
            'sample_size'])
        if 'differr' in tools:
            final_output += expand("results/differr/{comp}-{sample_size}/{comp}.differr.bed",comp=comparisons,sample_size=config[
            'sample_size'])
    else:
        if 'baleen' in tools:
            final_output += expand('results/baleen/{comp}/transcripts.csv',comp=comparisons)
        if 'nanocompore' in tools:
            if config['group']:
                final_output += [f"results/nanocompore/Group_{native_list}_{control_list}"]
            final_output += expand("results/modifications/{comp}/nanocompore.tsv.gz",comp=comparisons)
        if 'xpore' in tools:
            if 'genome' in config['params']['xpore']:
                final_output += expand("results/xpore/genome/{comp}/majority_direction_kmer_diffmod.table",comp=comparisons)
                if config['group']:
                    final_output += [f"results/xpore/Groups_genome/{native_list}_{control_list}/majority_direction_kmer_diffmod.table"]
            if config['group']:
                final_output += [f"results/xpore/Groups/{native_list}_{control_list}/majority_direction_kmer_diffmod.table"]
            final_output += expand("results/xpore/{comp}/majority_direction_kmer_diffmod.table",comp=comparisons)
        if 'm6anet' in tools:
            final_output += expand("results/m6anet/{sample}/data.site_proba.csv",sample=list(samples.index))
        if 'psinanopore' in tools:
            final_output += expand("results/psinanopore/{comp}.psi_candidates.csv",comp=comparisons)
        if 'differr' in tools:
            final_output += expand('results/differr/{comp}/{comp}.differr.bed',comp=comparisons)
    return final_output



