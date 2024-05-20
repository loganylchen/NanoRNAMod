import os
import sys
import pandas as pd

def format_xpore(input_file,output_file):
    df = pd.read_csv(input_file,sep='\t').sort_values(['id','position'])
    df.to_csv(output_file,sep='\t',index=False)

def format_baleen(input_file,output_file):
    df = pd.read_csv(input_file).sort_values(['transcript','loc'])
    df.to_csv(output_file,sep='\t',index=False)

def format_differr(input_file,output_file):
    columns = ['chrom', 'start', 'end', 'name', 'score',
               'strand', 'odds ratio',
               'G statistic ',
               '-log10 P value',
               '-log10 FDR',
               'G statistic for the homogeneity test of mutant replicates',
               '-log10 P value for the homogeneity test of mutant replicates',
               'G statistic for the homogeneity test of wild type replicates',
               '-log10 P value for the homogeneity test of wild type replicates', ]
    df = pd.read_csv(input_file, sep='\t', header=None)
    df.columns = columns
    df.to_csv(output_file,sep='\t',index=False)

if snakemake.params.tool == 'xpore':
    format_xpore(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'baleen':
    format_baleen(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'differr':
    format_differr(snakemake.input[0],snakemake.output[0])