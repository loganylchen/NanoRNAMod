import os
import sys
import pandas as pd
import glob

def format_xpore(input_file,output_file):
    df = pd.read_csv(input_file).sort_values(['id','position'])
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

def format_eligos2(input_file,output_file):
    df = pd.read_csv(input_file,sep='\t').sort_values(['chrom','start_loc'])
    df.to_csv(output_file,sep='\t',index=False)

def format_epinano(input_file,output_file):
    df = pd.read_csv(input_file)
    df[['chrom','pos','ref','strand']] = df['chr_pos'].str.split(' ',expand=True)
    df = df.loc[:,['chrom','pos','ref','strand','ko_feature',  'wt_feature' , 'delta_sum_err' , 'z_scores', 'z_score_prediction']].sort_values(['chrom','pos'])
    df.to_csv(output_file,sep='\t',index=False)

def format_drummer(input_dir,output_file):
    control = snakemake.wildcards.control
    native = snakemake.wildcards.native
    f = glob.glob(f'{input_dir}/{control}*-{native}*/summary.txt')[0]
    df = pd.read_csv(f,sep='\t').sort_values(['transcript_id','transcript_pos'])
    df.to_csv(output_file,sep='\t',index=False)


if snakemake.params.tool == 'xpore':
    format_xpore(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'baleen':
    format_baleen(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'differr':
    format_differr(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'epinano':
    format_epinano(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'eligos2':
    format_eligos2(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'drummer':
    format_drummer(snakemake.input[0],snakemake.output[0])