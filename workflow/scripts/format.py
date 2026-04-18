import os
import sys
import pandas as pd
import glob


def _input_has_content(path):
    """Return True if path points to a non-empty file.

    Used to short-circuit each format_* function when the upstream tool
    ran successfully but produced no results (empty-file safety net).
    """
    return os.path.isfile(path) and os.path.getsize(path) > 0


def _write_empty(output_file):
    """Create a zero-byte output file signalling 'no results'."""
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    open(output_file, "w").close()


def format_xpore(input_file,output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file).sort_values(['id','position'])
    df.to_csv(output_file,sep='\t',index=False)

def format_baleen(input_file,output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file).sort_values(['transcript','position'])
    df.to_csv(output_file,sep='\t',index=False)

def format_differr(input_file,output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
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
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file,sep='\t').sort_values(['chrom','start_loc'])
    df.to_csv(output_file,sep='\t',index=False)

def format_epinano(input_file,output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file)
    df[['chrom','pos','ref','strand']] = df['chr_pos'].str.split(' ',expand=True)
    # convert pos to 0-based
    df['pos'] = df['pos'].astype(int)-1
    df = df.loc[:,['chrom','pos','ref','strand','ko_feature',  'wt_feature' , 'delta_sum_err' , 'z_scores', 'z_score_prediction']].sort_values(['chrom','pos'])

    df.to_csv(output_file,sep='\t',index=False)

def format_drummer(input_dir,output_file):
    control = snakemake.wildcards.control
    native = snakemake.wildcards.native
    files = glob.glob(f'{input_dir}/{control}*-{native}*/summary.txt')
    # Tool succeeded but produced no summary (e.g. all sites filtered out by
    # depth/quality) → emit empty sentinel instead of failing.
    if not files:
        _write_empty(output_file); return
    f = files[0]
    if not _input_has_content(f):
        _write_empty(output_file); return
    df = pd.read_csv(f,sep='\t').sort_values(['transcript_id','transcript_pos'])
    # convert pos to 0-based
    df['transcript_pos'] = df['transcript_pos']-1
    df.to_csv(output_file,sep='\t',index=False)

def format_nanocompore(input_file,output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file,sep='\t').sort_values(['ref_id','pos'])
    # move the position by 2
    df['pos']+=2
    df.to_csv(output_file,sep='\t',index=False)

def format_tandemmod(input_file,output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file, sep='\t')
    if 'transcript' in df.columns and 'position' in df.columns:
        df = df.sort_values(['transcript', 'position'])
    elif 'id' in df.columns and 'position' in df.columns:
        df = df.sort_values(['id', 'position'])
        df.rename(columns={'id': 'transcript'}, inplace=True)
    df.to_csv(output_file, sep='\t', index=False)

def format_directrm(input_file,output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file, sep='\t')
    if 'transcript' in df.columns and 'position' in df.columns:
        df = df.sort_values(['transcript', 'position'])
    elif 'transcript_id' in df.columns and 'transcript_pos' in df.columns:
        df = df.sort_values(['transcript_id', 'transcript_pos'])
        df.rename(columns={'transcript_id': 'transcript', 'transcript_pos': 'position'}, inplace=True)
    df.to_csv(output_file, sep='\t', index=False)

def format_m6atm(input_file,output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file, sep='\t')
    if 'transcript' in df.columns and 'position' in df.columns:
        df = df.sort_values(['transcript', 'position'])
    elif 'transcript_id' in df.columns and 'transcript_pos' in df.columns:
        df = df.sort_values(['transcript_id', 'transcript_pos'])
        df.rename(columns={'transcript_id': 'transcript', 'transcript_pos': 'position'}, inplace=True)
    df.to_csv(output_file, sep='\t', index=False)

def format_rnano(input_file,output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file, sep='\t')
    if 'transcript' in df.columns and 'position' in df.columns:
        df = df.sort_values(['transcript', 'position'])
    elif 'transcript_id' in df.columns and 'transcript_pos' in df.columns:
        df = df.sort_values(['transcript_id', 'transcript_pos'])
        df.rename(columns={'transcript_id': 'transcript', 'transcript_pos': 'position'}, inplace=True)
    df.to_csv(output_file, sep='\t', index=False)

def format_psipore(input_file, output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file, sep='\t')
    if 'transcript' in df.columns and 'position' in df.columns:
        df = df.sort_values(['transcript', 'position'])
    df.to_csv(output_file, sep='\t', index=False)

def format_nanopsu(input_file, output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file, sep='\t')
    if 'transcript' in df.columns and 'position' in df.columns:
        df = df.sort_values(['transcript', 'position'])
    df.to_csv(output_file, sep='\t', index=False)

def format_nanomud(input_file, output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file, sep='\t')
    if 'transcript' in df.columns and 'position' in df.columns:
        sort_cols = ['transcript', 'position']
        if 'modification' in df.columns:
            sort_cols.append('modification')
        df = df.sort_values(sort_cols)
    df.to_csv(output_file, sep='\t', index=False)

def format_penguin(input_file, output_file):
    if not _input_has_content(input_file):
        _write_empty(output_file); return
    df = pd.read_csv(input_file, sep='\t')
    if 'transcript' in df.columns and 'position' in df.columns:
        df = df.sort_values(['transcript', 'position'])
    df.to_csv(output_file, sep='\t', index=False)

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
elif snakemake.params.tool == 'nanocompore':
    format_nanocompore(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'tandemmod':
    format_tandemmod(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'directrm':
    format_directrm(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'm6atm':
    format_m6atm(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'rnano':
    format_rnano(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'psipore':
    format_psipore(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'nanopsu':
    format_nanopsu(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'nanomud':
    format_nanomud(snakemake.input[0],snakemake.output[0])
elif snakemake.params.tool == 'penguin':
    format_penguin(snakemake.input[0],snakemake.output[0])
