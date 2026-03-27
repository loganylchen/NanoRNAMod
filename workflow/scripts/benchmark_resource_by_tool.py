"""
Aggregate and visualize FULL end-to-end resource usage by modification tool.

Each tool gets the COMPLETE cost of all its prerequisites (alignment, eventalign,
dataprep, etc.) plus its own execution. Shared prerequisites (minimap2, eventalign)
are counted in every tool that uses them, but tagged as "shared" so visualizations
can distinguish shared vs unique costs.

Tool dependency map (benchmark file stems):
- minimap2_genome_alignment, minimap2_transcriptome_alignment: shared across all tools
- f5c_index: shared across all eventalign-dependent tools + pybaleen
- eventalign_full: shared across xpore, nanocompore, baleen
- eventalign_simple: shared across psipore, nanopsu, nanomud
- Tool-specific prep/execution steps: unique to each tool
"""

import os
import glob
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
except ImportError:
    print("Error: matplotlib/seaborn not available")
    raise


# Metrics columns from Snakemake benchmark files
METRICS = ['s', 'cpu_time', 'max_rss', 'max_vms', 'io_in', 'io_out', 'mean_load']

# For max_rss and max_vms, we take max instead of sum across steps
MAX_METRICS = {'max_rss', 'max_vms'}

# Complete prerequisite chains for every tool.
# Each entry is (benchmark_stem, step_type) where step_type is 'shared' or 'unique'.
# benchmark_stem is the part between {sample}. and .benchmark.txt (or .txt).
#
# "shared" = used by multiple tools (alignment, eventalign)
# "unique" = specific to this tool's pipeline
#
# Note: link_fastq/link_blow5 are instant symlinks with no benchmark files.

TOOL_PREREQ_CHAINS = {
    # --- Comparison tools (run on {native}_{control}) ---
    'xpore': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('f5c_index', 'shared'),
        ('eventalign_full', 'shared'),
        ('uncompress_eventalign', 'unique'),
        ('xpore_dataprep', 'unique'),
        ('xpore', 'unique'),
        ('xpore_postprocessing', 'unique'),
    ],
    'nanocompore': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('f5c_index', 'shared'),
        ('eventalign_full', 'shared'),
        ('full_uncompress_eventalign', 'unique'),
        ('nanocompore_collapse', 'unique'),
        ('nanocompore', 'unique'),
    ],
    'baleen': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('f5c_index', 'shared'),
        ('eventalign_full', 'shared'),
        ('baleen_dataprep', 'unique'),
        ('baleen_modcall', 'unique'),
        ('baleen_postcall', 'unique'),
    ],
    'pybaleen': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('f5c_index', 'shared'),
        ('pybaleen', 'unique'),
    ],
    'differr': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('differr', 'unique'),
    ],
    'drummer': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('drummer', 'unique'),
    ],
    'eligos2': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('eligos2', 'unique'),
    ],
    'epinano': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('epinano_prep', 'unique'),
        ('epinano', 'unique'),
    ],
    'psipore': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('f5c_index', 'shared'),
        ('eventalign_simple', 'shared'),
        ('psipore', 'unique'),
    ],
    # --- Comparison tools that are also per-sample for eventalign ---
    'nanopsu': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('f5c_index', 'shared'),
        ('eventalign_simple', 'shared'),
        ('nanopsu', 'unique'),
    ],
    'nanomud': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('f5c_index', 'shared'),
        ('eventalign_simple', 'shared'),
        ('nanomud', 'unique'),
    ],
    # --- Per-sample tools ---
    'tandemmod': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('tandemmod', 'unique'),
    ],
    'directrm': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('directrm', 'unique'),
    ],
    'm6atm': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('m6atm', 'unique'),
    ],
    'rnano': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('rnano', 'unique'),
    ],
    'penguin': [
        ('minimap2_genome_alignment', 'shared'),
        ('minimap2_transcriptome_alignment', 'shared'),
        ('penguin', 'unique'),
    ],
}

# Tools where the main rule runs per-comparison ({native}_{control})
COMPARISON_TOOLS = {
    'xpore', 'nanocompore', 'baleen', 'pybaleen', 'differr', 'drummer',
    'eligos2', 'epinano', 'psipore',
}

# Tools where the main rule runs per-sample
PER_SAMPLE_TOOLS = {
    'tandemmod', 'directrm', 'm6atm', 'rnano', 'penguin', 'nanopsu', 'nanomud',
}

# Steps that always run per-sample regardless of tool type
PER_SAMPLE_STEPS = {
    'minimap2_genome_alignment', 'minimap2_transcriptome_alignment',
    'minimap2_transcriptome_alignment_3.2.4',
    'eventalign_full', 'eventalign_simple', 'f5c_index',
    'uncompress_eventalign', 'full_uncompress_eventalign',
    'xpore_dataprep', 'nanocompore_collapse', 'baleen_dataprep',
    'epinano_prep',
}


def find_benchmark_files(benchmark_dir):
    """Find all benchmark files and return a dict mapping (stem, sample_or_comparison) to file path."""
    files = {}
    for pattern in ['*.benchmark.txt', '*.txt']:
        for f in glob.glob(os.path.join(benchmark_dir, pattern)):
            basename = os.path.basename(f)
            # Extract stem: everything after first dot, minus the suffix
            if basename.endswith('.benchmark.txt'):
                without_suffix = basename[:-len('.benchmark.txt')]
            elif basename.endswith('.txt'):
                without_suffix = basename[:-len('.txt')]
            else:
                continue

            parts = without_suffix.rsplit('.', 1)
            if len(parts) != 2:
                continue

            sample_or_comp, stem = parts
            # Avoid duplicates (.benchmark.txt files also match *.txt)
            key = (stem, sample_or_comp)
            if key not in files:
                files[key] = f

    return files


def read_benchmark_file(filepath):
    """Read a single benchmark file and return metrics as a dict."""
    try:
        df = pd.read_csv(filepath, sep='\t')
        if df.empty:
            return None
        # Take the last row (Snakemake may write multiple if restarted)
        row = df.iloc[-1]
        return {m: row[m] if m in row.index else 0.0 for m in METRICS}
    except Exception as e:
        print(f"Warning: could not parse {filepath}: {e}")
        return None


def aggregate_by_tool(benchmark_dir):
    """
    For each tool, find all matching benchmark files for each prerequisite step,
    and build per-step breakdown + total.

    Returns a list of dicts, each representing one row in the output TSV.
    """
    file_index = find_benchmark_files(benchmark_dir)

    if not file_index:
        print("No benchmark files found")
        return []

    # Report what stems we found
    stems_found = sorted(set(stem for stem, _ in file_index.keys()))
    print(f"Benchmark stems found: {stems_found}")

    # Get all sample/comparison identifiers
    all_identifiers = sorted(set(sc for _, sc in file_index.keys()))
    print(f"Samples/comparisons found: {all_identifiers}")

    rows = []

    for tool, chain in TOOL_PREREQ_CHAINS.items():
        # Determine which identifiers this tool applies to
        # The tool's own benchmark stem tells us if it's per-sample or per-comparison
        tool_stem = tool  # The tool's own execution step has the same name as the tool
        tool_identifiers = [sc for (stem, sc) in file_index.keys() if stem == tool_stem]

        if not tool_identifiers:
            # Try alternate stems for tools with different benchmark names
            # (e.g., baleen uses baleen_modcall, not baleen)
            alt_stems = [stem for stem, stype in chain if stype == 'unique']
            for alt in alt_stems:
                tool_identifiers = [sc for (stem, sc) in file_index.keys() if stem == alt]
                if tool_identifiers:
                    break

        if not tool_identifiers:
            print(f"No benchmark data found for tool: {tool}")
            continue

        for identifier in tool_identifiers:
            # For comparison tools, identifier is like "A_B"
            # Per-sample steps need individual sample names
            if '_' in identifier and tool in COMPARISON_TOOLS:
                # Split comparison into native and control samples
                parts = identifier.split('_', 1)
                samples_for_shared = parts  # Both samples contribute to shared steps
                is_comparison = True
            else:
                samples_for_shared = [identifier]
                is_comparison = False

            for step_stem, step_type in chain:
                is_per_sample_step = step_stem in PER_SAMPLE_STEPS

                if is_per_sample_step and is_comparison:
                    # Sum across both samples in the comparison
                    step_metrics = {m: 0.0 for m in METRICS}
                    found_any = False
                    for sample in samples_for_shared:
                        key = (step_stem, sample)
                        if key in file_index:
                            m = read_benchmark_file(file_index[key])
                            if m:
                                found_any = True
                                for metric in METRICS:
                                    if metric in MAX_METRICS:
                                        step_metrics[metric] = max(step_metrics[metric], m[metric])
                                    else:
                                        step_metrics[metric] += m[metric]
                    if not found_any:
                        continue
                else:
                    # Look up directly with the identifier
                    key = (step_stem, identifier)
                    if key not in file_index:
                        continue
                    m = read_benchmark_file(file_index[key])
                    if not m:
                        continue
                    step_metrics = m

                row = {
                    'tool': tool,
                    'step': step_stem,
                    'step_type': step_type,
                    'sample_or_comparison': identifier,
                }
                row.update(step_metrics)
                rows.append(row)

    return rows


def build_output_dataframe(rows):
    """
    Build the output dataframe with per-step rows and summary total rows per tool.
    """
    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)

    # Add summary rows (one per tool per sample/comparison)
    summary_rows = []
    for (tool, identifier), group in df.groupby(['tool', 'sample_or_comparison']):
        summary = {
            'tool': tool,
            'step': 'TOTAL',
            'step_type': 'total',
            'sample_or_comparison': identifier,
        }
        for metric in METRICS:
            if metric in MAX_METRICS:
                summary[metric] = group[metric].max()
            else:
                summary[metric] = group[metric].sum()
        summary_rows.append(summary)

    summary_df = pd.DataFrame(summary_rows)
    result = pd.concat([df, summary_df], ignore_index=True)

    # Sort: by tool, then sample, then step (with TOTAL last)
    def sort_key(row):
        step_order = 0 if row['step'] != 'TOTAL' else 1
        return (row['tool'], row['sample_or_comparison'], step_order, row['step'])

    result['_sort'] = result.apply(sort_key, axis=1)
    result = result.sort_values('_sort').drop('_sort', axis=1).reset_index(drop=True)

    return result


def create_resource_visualization(df, output_path):
    """Create PDF visualization with shared/unique cost breakdown."""
    if df.empty:
        print("No data to visualize")
        return

    # Work only with TOTAL rows for overview charts
    totals = df[df['step'] == 'TOTAL'].copy()

    if totals.empty:
        print("No total rows to visualize")
        return

    with PdfPages(output_path) as pdf:
        # --- Page 1: Full-chain cost overview (4-panel) ---
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('Full End-to-End Resource Usage by Tool\n(includes all prerequisite steps)',
                     fontsize=14, fontweight='bold')

        tool_order = totals.groupby('tool')['s'].median().sort_values(ascending=False).index

        # Wall time
        ax = axes[0, 0]
        plot_df = totals.copy()
        plot_df['s_min'] = plot_df['s'] / 60
        sns.boxplot(data=plot_df, x='tool', y='s_min', order=tool_order, ax=ax, palette='Set2')
        ax.set_xlabel('')
        ax.set_ylabel('Total Time (minutes)')
        ax.set_title('Wall Time (full chain)')
        ax.tick_params(axis='x', rotation=45)

        # CPU time
        ax = axes[0, 1]
        plot_df['cpu_min'] = plot_df['cpu_time'] / 60
        sns.boxplot(data=plot_df, x='tool', y='cpu_min', order=tool_order, ax=ax, palette='Set2')
        ax.set_xlabel('')
        ax.set_ylabel('CPU Time (minutes)')
        ax.set_title('CPU Time (full chain)')
        ax.tick_params(axis='x', rotation=45)

        # Memory
        ax = axes[1, 0]
        plot_df['rss_gb'] = plot_df['max_rss'] / 1024
        sns.boxplot(data=plot_df, x='tool', y='rss_gb', order=tool_order, ax=ax, palette='Set2')
        ax.set_xlabel('')
        ax.set_ylabel('Peak RSS (GB)')
        ax.set_title('Peak Memory')
        ax.tick_params(axis='x', rotation=45)

        # I/O
        ax = axes[1, 1]
        plot_df['io_gb'] = (plot_df['io_in'] + plot_df['io_out']) / 1024
        sns.boxplot(data=plot_df, x='tool', y='io_gb', order=tool_order, ax=ax, palette='Set2')
        ax.set_xlabel('')
        ax.set_ylabel('I/O (GB)')
        ax.set_title('Total I/O')
        ax.tick_params(axis='x', rotation=45)

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # --- Page 2: Stacked bar with shared/unique breakdown ---
        fig, ax = plt.subplots(figsize=(14, 8))

        # Per-step rows (not TOTAL)
        steps_df = df[df['step'] != 'TOTAL'].copy()

        # Aggregate: mean wall time per tool per step_type
        breakdown = steps_df.groupby(['tool', 'step_type'])['s'].mean().unstack(fill_value=0) / 60
        # Ensure both columns exist
        for col in ['shared', 'unique']:
            if col not in breakdown.columns:
                breakdown[col] = 0.0
        breakdown = breakdown[['shared', 'unique']]
        breakdown = breakdown.loc[tool_order]

        colors = {'shared': '#95a5a6', 'unique': '#3498db'}
        breakdown.plot(kind='bar', stacked=True, ax=ax,
                       color=[colors['shared'], colors['unique']])

        ax.set_xlabel('Tool')
        ax.set_ylabel('Mean Time (minutes)')
        ax.set_title('Time Breakdown: Shared vs Unique Steps\n'
                      '(Shared = alignment, eventalign, f5c_index; Unique = tool-specific)',
                      fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.legend(title='Step Type', labels=['Shared prerequisites', 'Unique to tool'])

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # --- Page 3: Detailed step breakdown (top tools) ---
        fig, ax = plt.subplots(figsize=(16, 9))

        # Pivot: mean time per tool per step
        step_times = steps_df.groupby(['tool', 'step'])['s'].mean().unstack(fill_value=0) / 60

        # Get all unique steps in chain order
        all_steps = []
        for chain in TOOL_PREREQ_CHAINS.values():
            for stem, _ in chain:
                if stem not in all_steps:
                    all_steps.append(stem)

        # Only include steps that have data
        available_steps = [s for s in all_steps if s in step_times.columns]
        step_times = step_times.reindex(columns=available_steps, fill_value=0)
        step_times = step_times.loc[step_times.index.isin(tool_order)]
        step_times = step_times.loc[tool_order]

        # Color: shared steps in gray tones, unique in color
        shared_steps = set()
        for chain in TOOL_PREREQ_CHAINS.values():
            for stem, stype in chain:
                if stype == 'shared':
                    shared_steps.add(stem)

        cmap_unique = plt.cm.tab20
        cmap_shared = plt.cm.Greys
        step_colors = []
        unique_idx = 0
        shared_idx = 0
        for step in available_steps:
            if step in shared_steps:
                step_colors.append(cmap_shared(0.3 + shared_idx * 0.15))
                shared_idx += 1
            else:
                step_colors.append(cmap_unique(unique_idx * 0.06 + 0.1))
                unique_idx += 1

        step_times.plot(kind='bar', stacked=True, ax=ax, color=step_colors)
        ax.set_xlabel('Tool')
        ax.set_ylabel('Mean Time (minutes)')
        ax.set_title('Detailed Step Breakdown per Tool\n(gray = shared steps, color = unique steps)',
                      fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.legend(title='Step', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7)

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # --- Page 4: Time vs Memory scatter ---
        fig, ax = plt.subplots(figsize=(12, 8))

        summary = totals.groupby('tool').agg({
            's': 'mean',
            'max_rss': 'mean',
        }).reset_index()

        summary['s_min'] = summary['s'] / 60
        summary['rss_gb'] = summary['max_rss'] / 1024

        # Color by number of shared steps
        def count_shared(tool):
            chain = TOOL_PREREQ_CHAINS.get(tool, [])
            return sum(1 for _, t in chain if t == 'shared')

        summary['n_shared'] = summary['tool'].apply(count_shared)

        scatter = ax.scatter(summary['s_min'], summary['rss_gb'],
                            c=summary['n_shared'], cmap='RdYlGn_r',
                            s=200, alpha=0.8, edgecolors='black', vmin=2, vmax=5)

        for _, row in summary.iterrows():
            ax.annotate(row['tool'], (row['s_min'], row['rss_gb']),
                       xytext=(5, 5), textcoords='offset points', fontsize=9)

        plt.colorbar(scatter, ax=ax, label='Number of shared prerequisite steps')
        ax.set_xlabel('Mean Total Wall Time (minutes)')
        ax.set_ylabel('Mean Peak Memory (GB)')
        ax.set_title('Resource Profile: Full-Chain Time vs Peak Memory',
                     fontweight='bold')
        ax.grid(alpha=0.3)

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # --- Page 5: Summary table ---
        fig, ax = plt.subplots(figsize=(16, 10))
        ax.axis('off')

        table_data = totals.groupby('tool').agg({
            's': ['mean', 'std'],
            'cpu_time': ['mean', 'std'],
            'max_rss': ['mean', 'max'],
            'io_in': 'mean',
            'io_out': 'mean',
        }).round(1)

        table_data.columns = [
            'Time Mean (s)', 'Time Std (s)',
            'CPU Mean (s)', 'CPU Std (s)',
            'Mem Mean (MB)', 'Mem Max (MB)',
            'I/O In Mean (MB)', 'I/O Out Mean (MB)',
        ]
        table_data = table_data.reset_index().rename(columns={'tool': 'Tool'})

        table = ax.table(cellText=table_data.values,
                        colLabels=table_data.columns,
                        cellLoc='center',
                        loc='center')

        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.2, 1.5)

        for i in range(len(table_data.columns)):
            table[(0, i)].set_facecolor('#2c3e50')
            table[(0, i)].set_text_props(color='white', fontweight='bold')

        for i in range(1, len(table_data) + 1):
            for j in range(len(table_data.columns)):
                if i % 2 == 0:
                    table[(i, j)].set_facecolor('#ecf0f1')

        ax.set_title('Full-Chain Resource Summary (all prerequisites included)',
                     fontsize=14, fontweight='bold', y=0.98)

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    print(f"Resource visualization saved to {output_path}")


def main():
    if 'snakemake' in globals():
        benchmark_dir = snakemake.input.benchmark_dir if hasattr(snakemake.input, 'benchmark_dir') else snakemake.params.benchmark_dir
        output_tsv = snakemake.output.tsv
        output_pdf = snakemake.output.pdf
    else:
        parser = argparse.ArgumentParser(description='Aggregate full-chain resources by tool')
        parser.add_argument('--benchmark-dir', required=True, help='Directory with benchmark files')
        parser.add_argument('--output-tsv', required=True, help='Output TSV file')
        parser.add_argument('--output-pdf', required=True, help='Output PDF visualization')
        args = parser.parse_args()

        benchmark_dir = args.benchmark_dir
        output_tsv = args.output_tsv
        output_pdf = args.output_pdf

    print(f"Loading benchmark files from {benchmark_dir}")
    rows = aggregate_by_tool(benchmark_dir)

    if not rows:
        print("No benchmark data found")
        pd.DataFrame().to_csv(output_tsv, sep='\t', index=False)
        return

    df = build_output_dataframe(rows)
    print(f"Generated {len(df)} rows for {df['tool'].nunique()} tools")

    # Save TSV
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Saved to {output_tsv}")

    # Create visualization
    print("Creating visualization...")
    create_resource_visualization(df, output_pdf)


if __name__ == '__main__':
    main()
