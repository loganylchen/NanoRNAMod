import os
import pandas as pd
import numpy as np


def normalize_columns(df):
    col_mapping = {}

    transcript_cols = ['transcript', 'id', 'ref_id', 'chrom', 'transcript_id']
    for col in transcript_cols:
        if col in df.columns:
            col_mapping[col] = 'transcript'
            break

    position_cols = ['position', 'pos', 'start', 'start_loc', 'transcript_pos',
                     'transcript_loc', 'genomic_pos']
    for col in position_cols:
        if col in df.columns and col not in col_mapping.values():
            col_mapping[col] = 'position'
            break

    if col_mapping:
        df = df.rename(columns=col_mapping)

    return df


def parse_windows(window_param):
    if isinstance(window_param, list):
        return [int(w) for w in window_param]
    if isinstance(window_param, (int, float)):
        return [int(window_param)]
    s = str(window_param).strip()
    if s.startswith('[') and s.endswith(']'):
        inner = s[1:-1]
        parts = [p.strip() for p in inner.split(',') if p.strip()]
        return [int(p) for p in parts]
    return [int(s)]


def run_per_tool_mode(truth_path, results_path, output_covered, window_param):
    windows = parse_windows(window_param)
    max_window = max(windows)

    truth = pd.read_csv(truth_path, sep='\t')
    if 'label' in truth.columns:
        truth = truth[truth['label'] != '-'].copy()
    truth = normalize_columns(truth)

    if truth.empty or 'transcript' not in truth.columns or 'position' not in truth.columns:
        pd.DataFrame(columns=['transcript', 'position']).to_csv(
            output_covered, sep='\t', index=False
        )
        return

    sep = ',' if results_path.endswith('.csv') else '\t'
    try:
        preds = pd.read_csv(results_path, sep=sep)
    except Exception:
        pd.DataFrame(columns=['transcript', 'position']).to_csv(
            output_covered, sep='\t', index=False
        )
        return

    preds = normalize_columns(preds)

    if preds.empty or 'transcript' not in preds.columns or 'position' not in preds.columns:
        pd.DataFrame(columns=['transcript', 'position']).to_csv(
            output_covered, sep='\t', index=False
        )
        return

    preds['position'] = pd.to_numeric(preds['position'], errors='coerce')
    preds = preds.dropna(subset=['position'])
    preds['position'] = preds['position'].astype(int)

    pred_by_tx = {}
    for tx, group in preds.groupby('transcript'):
        pred_by_tx[tx] = group['position'].values

    covered_rows = []
    for _, row in truth.iterrows():
        tx = row['transcript']
        pos = int(row['position'])
        if tx not in pred_by_tx:
            continue
        positions = pred_by_tx[tx]
        if np.any(np.abs(positions - pos) <= max_window):
            covered_rows.append({'transcript': tx, 'position': pos})

    out_df = pd.DataFrame(covered_rows, columns=['transcript', 'position'])
    out_df.to_csv(output_covered, sep='\t', index=False)


def run_union_mode(covered_paths, output_union, output_called_sites):
    site_tools = {}

    for path in covered_paths:
        filename = os.path.basename(path)
        if filename.endswith('_covered.tsv'):
            tool = filename[: -len('_covered.tsv')]
        else:
            tool = filename.replace('.tsv', '')

        try:
            df = pd.read_csv(path, sep='\t')
        except pd.errors.EmptyDataError:
            continue
        except Exception:
            continue

        if df.empty or 'transcript' not in df.columns or 'position' not in df.columns:
            continue

        df['position'] = pd.to_numeric(df['position'], errors='coerce')
        df = df.dropna(subset=['position'])
        df['position'] = df['position'].astype(int)

        for _, row in df.iterrows():
            key = (row['transcript'], int(row['position']))
            if key not in site_tools:
                site_tools[key] = set()
            site_tools[key].add(tool)

    if not site_tools:
        pd.DataFrame(columns=['transcript', 'position', 'covered_by_tools']).to_csv(
            output_union, sep='\t', index=False
        )
        pd.DataFrame(columns=['tool', 'n_covered']).to_csv(
            output_called_sites, sep='\t', index=False
        )
        return

    union_rows = []
    for (tx, pos), tools_set in sorted(site_tools.items()):
        union_rows.append({
            'transcript': tx,
            'position': pos,
            'covered_by_tools': ','.join(sorted(tools_set)),
        })

    union_df = pd.DataFrame(union_rows)
    union_df.to_csv(output_union, sep='\t', index=False)

    tool_counts = {}
    for tools_set in site_tools.values():
        for tool in tools_set:
            tool_counts[tool] = tool_counts.get(tool, 0) + 1

    called_rows = [{'tool': t, 'n_covered': n} for t, n in sorted(tool_counts.items())]
    called_df = pd.DataFrame(called_rows, columns=['tool', 'n_covered'])
    called_df.to_csv(output_called_sites, sep='\t', index=False)


output_keys = list(snakemake.output.keys())

if 'covered' in output_keys:
    run_per_tool_mode(
        truth_path=snakemake.input.truth_set,
        results_path=snakemake.input.results,
        output_covered=snakemake.output.covered,
        window_param=snakemake.params.window,
    )
elif 'union' in output_keys:
    covered_paths = list(snakemake.input.covered)
    run_union_mode(
        covered_paths=covered_paths,
        output_union=snakemake.output.union,
        output_called_sites=snakemake.output.called_sites,
    )
else:
    raise ValueError(
        "Cannot determine mode: snakemake.output must have 'covered' or 'union' key"
    )
