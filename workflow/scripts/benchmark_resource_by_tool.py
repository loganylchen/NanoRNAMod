"""
Aggregate and visualize resource usage by modification tool.

This script:
1. Maps benchmark files to modification tools
2. Includes prerequisite steps (alignment, eventalign) for each tool
3. Calculates total resources per tool
4. Generates visualizations showing resource usage over time

Tool dependency map:
- All tools need: link_fastq, minimap2 alignment
- Signal-based tools need: f5c eventalign (nanocompore, baleen, xpore, pybaleen, psipore, nanopsu, nanomud)
- Alignment-based tools need: just BAM (differr, drummer, eligos2, epinano)
- Per-sample tools: tandemmod, directrm, m6atm, rnano, penguin
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


# Define tool categories and their dependencies
COMPARISON_TOOLS = [
    'xpore', 'nanocompore', 'baleen', 'pybaleen', 'differr', 'drummer',
    'eligos2', 'epinano', 'psipore', 'nanopsu', 'nanomud'
]

PER_SAMPLE_TOOLS = ['tandemmod', 'directrm', 'm6atm', 'rnano', 'penguin']

# Tools that need eventalign (signal-based, use pre-processed eventalign TSV)
SIGNAL_TOOLS = ['nanocompore', 'baleen', 'xpore', 'psipore', 'nanopsu', 'nanomud']

# Tools that need f5c index (read raw signal files directly)
F5C_INDEX_TOOLS = ['pybaleen']

# Tools that only need alignment
ALIGNMENT_TOOLS = ['differr', 'drummer', 'eligos2', 'epinano', 'tandemmod', 'directrm', 'm6atm', 'rnano', 'penguin']

# Detailed prerequisites for each tool (for documentation page)
TOOL_PREREQUISITES = {
    'xpore': {
        'description': 'Differential RNA modification detection using nanopore direct RNA sequencing',
        'inputs': ['Native BAM (minimap2)', 'Control BAM (minimap2)', 'Eventalign TSV (f5c eventalign --collapse)', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2', 'f5c_eventalign'],
        'category': 'Signal-based (eventalign TSV)'
    },
    'nanocompore': {
        'description': 'Nanopore RNA modification calling using signal intensity and dwell time',
        'inputs': ['Native BAM (minimap2)', 'Control BAM (minimap2)', 'Eventalign TSV (f5c eventalign --collapse)', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2', 'f5c_eventalign'],
        'category': 'Signal-based (eventalign TSV)'
    },
    'baleen': {
        'description': 'RNA modification detection using k-mer based models',
        'inputs': ['Native BAM (minimap2)', 'Control BAM (minimap2)', 'Eventalign TSV (f5c eventalign --collapse)', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2', 'f5c_eventalign'],
        'category': 'Signal-based (eventalign TSV)'
    },
    'pybaleen': {
        'description': 'CUDA-accelerated DTW + HMM modification detection',
        'inputs': ['Native BAM', 'Control BAM', 'Native FASTQ', 'Control FASTQ', 'Native BLOW5 + Index', 'Control BLOW5 + Index', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'link_blow5', 'minimap2', 'f5c_index'],
        'category': 'Signal-based (raw BLOW5, needs GPU)'
    },
    'differr': {
        'description': 'Differential error rate analysis for RNA modification detection',
        'inputs': ['Native BAM (minimap2)', 'Control BAM (minimap2)', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2'],
        'category': 'Alignment-based'
    },
    'drummer': {
        'description': 'DRUMMER - Differential RNA Modification detection',
        'inputs': ['Native BAM (minimap2)', 'Control BAM (minimap2)', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2'],
        'category': 'Alignment-based'
    },
    'eligos2': {
        'description': 'Epitranscriptional landscape inferring from RNA sequencing',
        'inputs': ['Native BAM (minimap2)', 'Control BAM (minimap2)', 'Reference FASTA', 'Feature annotations (GTF)'],
        'prereq_rules': ['link_fastq', 'minimap2'],
        'category': 'Alignment-based'
    },
    'epinano': {
        'description': 'Detection of m6A RNA modifications using nanopore direct RNA sequencing',
        'inputs': ['Native BAM (minimap2)', 'Control BAM (minimap2)', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2'],
        'category': 'Alignment-based'
    },
    'psipore': {
        'description': 'RNA modification detection using PSI model',
        'inputs': ['Native BAM (minimap2)', 'Control BAM (minimap2)', 'Eventalign TSV (f5c eventalign --collapse)', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2', 'f5c_eventalign'],
        'category': 'Signal-based (eventalign TSV)'
    },
    'nanopsu': {
        'description': 'RNA modification detection using PSU features',
        'inputs': ['Native BAM (minimap2)', 'Control BAM (minimap2)', 'Eventalign TSV (f5c eventalign --collapse)', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2', 'f5c_eventalign'],
        'category': 'Signal-based (eventalign TSV)'
    },
    'nanomud': {
        'description': 'RNA modification detection using MUD model',
        'inputs': ['Native BAM (minimap2)', 'Control BAM (minimap2)', 'Eventalign TSV (f5c eventalign --collapse)', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2', 'f5c_eventalign'],
        'category': 'Signal-based (eventalign TSV)'
    },
    'tandemmod': {
        'description': 'RNA modification detection without control sample',
        'inputs': ['Native BAM (minimap2)', 'Native FASTQ', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2'],
        'category': 'Per-sample (no control needed)'
    },
    'directrm': {
        'description': 'Direct RNA modification detection',
        'inputs': ['Native BAM (minimap2)', 'Native FASTQ', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2'],
        'category': 'Per-sample (no control needed)'
    },
    'm6atm': {
        'description': 'm6A detection using transformer model',
        'inputs': ['Native BAM (minimap2)', 'Native FASTQ', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2'],
        'category': 'Per-sample (no control needed)'
    },
    'rnano': {
        'description': 'RNA modification detection using RNN',
        'inputs': ['Native BAM (minimap2)', 'Native FASTQ', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2'],
        'category': 'Per-sample (no control needed)'
    },
    'penguin': {
        'description': 'RNA modification detection',
        'inputs': ['Native BAM (minimap2)', 'Native FASTQ', 'Reference FASTA'],
        'prereq_rules': ['link_fastq', 'minimap2'],
        'category': 'Per-sample (no control needed)'
    }
}


def parse_benchmark_filename(filename):
    """
    Parse benchmark filename to extract sample/comparison and rule type.

    Examples:
        SampleA.minimap2_genome.benchmark.txt -> sample='SampleA', rule='minimap2_genome'
        SampleA_SampleB.xpore.benchmark.txt -> comparison='SampleA_SampleB', rule='xpore'
        SampleA.f5c_eventalign.benchmark.txt -> sample='SampleA', rule='f5c_eventalign'
    """
    stem = os.path.basename(filename).replace('.benchmark.txt', '')
    parts = stem.rsplit('.', 1)

    if len(parts) == 2:
        sample_or_comparison = parts[0]
        rule = parts[1]
    else:
        sample_or_comparison = stem
        rule = 'unknown'

    # Determine if it's a comparison (contains underscore and matches tool pattern)
    is_comparison = '_' in sample_or_comparison and rule in COMPARISON_TOOLS

    return {
        'sample_or_comparison': sample_or_comparison,
        'rule': rule,
        'is_comparison': is_comparison
    }


def categorize_rule(rule_name):
    """
    Categorize a rule into: modification_tool, alignment, eventalign, f5c_index, or other.
    """
    rule_lower = rule_name.lower()

    # Check if it's a modification tool
    for tool in COMPARISON_TOOLS + PER_SAMPLE_TOOLS:
        if tool in rule_lower:
            return ('modification_tool', tool)

    # Check if it's alignment
    if any(x in rule_lower for x in ['minimap2', 'alignment', 'align']):
        return ('alignment', 'minimap2')

    # Check if it's f5c index (must check before eventalign since both contain 'f5c')
    if 'f5c_index' in rule_lower or ('f5c' in rule_lower and 'index' in rule_lower):
        return ('f5c_index', 'f5c')

    # Check if it's eventalign
    if any(x in rule_lower for x in ['eventalign', 'f5c']):
        return ('eventalign', 'f5c')

    # Check if it's a linking step
    if 'link' in rule_lower:
        return ('link', 'link')

    return ('other', rule_name)


def get_tool_dependencies(tool_name):
    """
    Get the prerequisite rule categories for a modification tool.

    Returns list of categories: ['link', 'alignment', ...] + additional deps based on tool type
    """
    dependencies = ['link', 'alignment']

    tool_lower = tool_name.lower()

    if tool_lower in [t.lower() for t in SIGNAL_TOOLS]:
        # Tools that use pre-processed eventalign TSV
        dependencies.append('eventalign')
    elif tool_lower in [t.lower() for t in F5C_INDEX_TOOLS]:
        # Tools that read raw signal files directly (need f5c index)
        dependencies.append('f5c_index')

    return dependencies


def load_benchmark_files(benchmark_dir):
    """Load all benchmark files from directory."""
    benchmark_files = glob.glob(os.path.join(benchmark_dir, '*.benchmark.txt'))

    records = []
    for f in benchmark_files:
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty:
                continue

            parsed = parse_benchmark_filename(f)
            category, subcategory = categorize_rule(parsed['rule'])

            df['file'] = os.path.basename(f)
            df['sample_or_comparison'] = parsed['sample_or_comparison']
            df['rule'] = parsed['rule']
            df['rule_category'] = category
            df['rule_subcategory'] = subcategory
            df['is_comparison'] = parsed['is_comparison']

            records.append(df)
        except Exception as e:
            print(f"Warning: could not parse {f}: {e}")

    if not records:
        return pd.DataFrame()

    return pd.concat(records, ignore_index=True)


def aggregate_by_tool_with_dependencies(df):
    """
    Aggregate resource usage by modification tool, including dependencies.

    For each modification tool, sum up:
    - The tool's own resource usage
    - Alignment resource usage (prerequisite)
    - Eventalign resource usage (prerequisite for signal-based tools)
    """
    if df.empty:
        return pd.DataFrame()

    # Group by sample/comparison and rule category
    results = []

    # Get unique modification tools and their samples/comparisons
    mod_df = df[df['rule_category'] == 'modification_tool']

    if mod_df.empty:
        return pd.DataFrame()

    for tool in mod_df['rule_subcategory'].unique():
        tool_df = mod_df[mod_df['rule_subcategory'] == tool]
        dependencies = get_tool_dependencies(tool)

        for sample_comp in tool_df['sample_or_comparison'].unique():
            # Get data for this sample/comparison
            sample_df = df[df['sample_or_comparison'] == sample_comp]

            row = {
                'tool': tool,
                'sample_or_comparison': sample_comp,
            }

            # Aggregate metrics
            metrics = ['s', 'cpu_time', 'max_rss', 'max_vms', 'io_in', 'io_out', 'mean_load']

            # Tool's own resources
            tool_only = sample_df[sample_df['rule_subcategory'] == tool]
            for metric in metrics:
                if metric in tool_only.columns:
                    row[f'{metric}_tool'] = tool_only[metric].sum()

            # Alignment resources
            align_df = sample_df[sample_df['rule_category'] == 'alignment']
            for metric in metrics:
                if metric in align_df.columns:
                    row[f'{metric}_alignment'] = align_df[metric].sum()

            # Eventalign resources (if applicable)
            if 'eventalign' in dependencies:
                event_df = sample_df[sample_df['rule_category'] == 'eventalign']
                for metric in metrics:
                    if metric in event_df.columns:
                        row[f'{metric}_eventalign'] = event_df[metric].sum()

            # f5c_index resources (if applicable - for tools that read raw signal directly)
            if 'f5c_index' in dependencies:
                index_df = sample_df[sample_df['rule_category'] == 'f5c_index']
                for metric in metrics:
                    if metric in index_df.columns:
                        row[f'{metric}_f5c_index'] = index_df[metric].sum()

            # Calculate total (tool + dependencies)
            for metric in metrics:
                total = 0
                for suffix in ['_tool', '_alignment', '_eventalign', '_f5c_index']:
                    key = f'{metric}{suffix}'
                    if key in row and pd.notna(row[key]):
                        total += row[key]
                row[f'{metric}_total'] = total

            results.append(row)

    return pd.DataFrame(results)


def create_prerequisites_page(pdf, agg_df):
    """Create a page showing prerequisites for each tool."""
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111)
    ax.axis('off')

    ax.text(0.5, 0.97, 'Tool Prerequisites Overview',
            fontsize=20, fontweight='bold', ha='center', va='center',
            transform=ax.transAxes)

    ax.text(0.5, 0.93, 'Input requirements and prerequisite steps for each modification detection tool',
            fontsize=11, ha='center', va='center', transform=ax.transAxes,
            color='gray')

    # Get tools from data
    tools_in_data = sorted(agg_df['tool'].unique()) if not agg_df.empty else []

    # Build table data
    table_data = []
    headers = ['Tool', 'Category', 'Prerequisite Steps', 'Key Inputs']

    for tool in tools_in_data:
        info = TOOL_PREREQUISITES.get(tool, {})
        category = info.get('category', 'Unknown')
        prereqs = info.get('prereq_rules', [])
        inputs = info.get('inputs', [])

        # Format for display
        prereqs_str = ', '.join(prereqs[:3])
        if len(prereqs) > 3:
            prereqs_str += f' +{len(prereqs)-3}'

        inputs_str = '\n'.join(inputs[:2])
        if len(inputs) > 2:
            inputs_str += f'\n+{len(inputs)-2} more'

        table_data.append([tool, category, prereqs_str, inputs_str])

    if not table_data:
        ax.text(0.5, 0.5, 'No tool data available',
                fontsize=12, ha='center', va='center', transform=ax.transAxes)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        return

    # Create table
    table = ax.table(cellText=table_data,
                     colLabels=headers,
                     cellLoc='left',
                     loc='center',
                     colWidths=[0.12, 0.20, 0.28, 0.40])

    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.1, 2.2)

    # Style header
    for i in range(len(headers)):
        table[(0, i)].set_facecolor('#2c3e50')
        table[(0, i)].set_text_props(color='white', fontweight='bold')

    # Alternate row colors
    for i in range(1, len(table_data) + 1):
        for j in range(len(headers)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#ecf0f1')

    # Set column alignments
    for i in range(len(table_data) + 1):
        for j in range(len(headers)):
            cell = table[(i, j)]
            cell.get_text().set_wrap(True)

    # Add legend for prerequisite steps
    legend_text = """Prerequisite Step Definitions:
• link_fastq: Symlink raw FASTQ files to project results
• link_blow5: Symlink raw BLOW5 signal files to project results
• minimap2: Align reads to reference transcriptome (produces BAM)
• f5c_index: Create BLOW5 index files for fast signal access
• f5c_eventalign: Align signal events to reference (produces TSV)

Category Legend:
• Signal-based: Requires processed signal data (eventalign TSV or raw BLOW5)
• Alignment-based: Only requires BAM alignment files
• Per-sample: No control sample needed"""

    ax.text(0.02, 0.15, legend_text,
            fontsize=7, va='bottom', ha='left',
            transform=ax.transAxes, family='monospace',
            bbox=dict(boxstyle='round', facecolor='#f8f9fa', alpha=0.8))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def create_resource_visualization(agg_df, output_path):
    """Create PDF visualization of resource usage by tool."""
    if agg_df.empty:
        print("No data to visualize")
        return

    with PdfPages(output_path) as pdf:
        # Page 1: Tool Prerequisites Overview
        create_prerequisites_page(pdf, agg_df)

        # Page 2: Total wall time by tool
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('Resource Usage by Modification Tool (Total including prerequisites)',
                     fontsize=14, fontweight='bold')

        # Sort tools by median total time
        tool_order = agg_df.groupby('tool')['s_total'].median().sort_values(ascending=False).index

        # 1. Wall time boxplot
        ax = axes[0, 0]
        plot_df = agg_df.copy()
        plot_df['s_total_min'] = plot_df['s_total'] / 60  # Convert to minutes
        sns.boxplot(data=plot_df, x='tool', y='s_total_min', order=tool_order, ax=ax, palette='Set2')
        ax.set_xlabel('Tool')
        ax.set_ylabel('Total Time (minutes)')
        ax.set_title('Total Wall Time (tool + prerequisites)')
        ax.tick_params(axis='x', rotation=45)

        # 2. CPU time boxplot
        ax = axes[0, 1]
        plot_df['cpu_time_min'] = plot_df['cpu_time_total'] / 60
        sns.boxplot(data=plot_df, x='tool', y='cpu_time_min', order=tool_order, ax=ax, palette='Set2')
        ax.set_xlabel('Tool')
        ax.set_ylabel('CPU Time (minutes)')
        ax.set_title('Total CPU Time')
        ax.tick_params(axis='x', rotation=45)

        # 3. Memory usage
        ax = axes[1, 0]
        plot_df['max_rss_gb'] = plot_df['max_rss_total'] / 1024  # Convert MB to GB
        sns.boxplot(data=plot_df, x='tool', y='max_rss_gb', order=tool_order, ax=ax, palette='Set2')
        ax.set_xlabel('Tool')
        ax.set_ylabel('Max RSS (GB)')
        ax.set_title('Peak Memory Usage')
        ax.tick_params(axis='x', rotation=45)

        # 4. I/O
        ax = axes[1, 1]
        plot_df['io_total_gb'] = (plot_df['io_in_total'] + plot_df['io_out_total']) / 1024
        sns.boxplot(data=plot_df, x='tool', y='io_total_gb', order=tool_order, ax=ax, palette='Set2')
        ax.set_xlabel('Tool')
        ax.set_ylabel('I/O (GB)')
        ax.set_title('Total I/O (read + write)')
        ax.tick_params(axis='x', rotation=45)

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 2: Time breakdown by component
        fig, ax = plt.subplots(figsize=(14, 8))

        # Calculate mean time by tool and component
        # Note: f5c_index is for tools that read raw signal directly (pybaleen)
        # eventalign is for tools that use pre-processed signal TSV
        components = ['tool', 'alignment', 'eventalign', 'f5c_index']
        time_cols = [f's_{c}' for c in components]

        # Filter to only columns that exist
        available_cols = [c for c in time_cols if c in agg_df.columns]
        mean_times = agg_df.groupby('tool')[available_cols].mean() / 60  # to minutes
        mean_times = mean_times.loc[tool_order]

        # Rename columns for display
        col_rename = {
            's_tool': 'Tool\nItself',
            's_alignment': 'Alignment\n(prereq)',
            's_eventalign': 'Eventalign\n(prereq)',
            's_f5c_index': 'F5C Index\n(prereq)'
        }
        mean_times = mean_times.rename(columns={k: v for k, v in col_rename.items() if k in mean_times.columns})

        # Colors for each component
        color_map = {
            'Tool\nItself': '#3498db',
            'Alignment\n(prereq)': '#95a5a6',
            'Eventalign\n(prereq)': '#e74c3c',
            'F5C Index\n(prereq)': '#f39c12'
        }
        colors = [color_map[c] for c in mean_times.columns if c in color_map]

        mean_times.plot(kind='bar', stacked=True, ax=ax, color=colors)
        ax.set_xlabel('Tool')
        ax.set_ylabel('Time (minutes)')
        ax.set_title('Time Breakdown: Tool vs Prerequisites', fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.legend(title='Component', loc='upper right')

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 3: Resource efficiency (time vs memory)
        fig, ax = plt.subplots(figsize=(12, 8))

        # Calculate mean per tool
        summary = agg_df.groupby('tool').agg({
            's_total': 'mean',
            'max_rss_total': 'mean',
            'cpu_time_total': 'mean'
        }).reset_index()

        summary['s_total_min'] = summary['s_total'] / 60
        summary['max_rss_gb'] = summary['max_rss_total'] / 1024

        # Color by whether tool needs signal processing (eventalign or f5c_index)
        def get_tool_color(tool_name):
            tool_lower = tool_name.lower()
            if tool_lower in [t.lower() for t in SIGNAL_TOOLS]:
                return '#e74c3c'  # Red for eventalign-based tools
            elif tool_lower in [t.lower() for t in F5C_INDEX_TOOLS]:
                return '#f39c12'  # Orange for f5c_index-based tools
            else:
                return '#3498db'  # Blue for alignment-only tools

        colors = [get_tool_color(t) for t in summary['tool']]

        scatter = ax.scatter(summary['s_total_min'], summary['max_rss_gb'],
                            c=colors, s=200, alpha=0.7, edgecolors='black')

        for i, row in summary.iterrows():
            ax.annotate(row['tool'], (row['s_total_min'], row['max_rss_gb']),
                       xytext=(5, 5), textcoords='offset points', fontsize=9)

        ax.set_xlabel('Mean Total Wall Time (minutes)')
        ax.set_ylabel('Mean Peak Memory (GB)')
        ax.set_title('Resource Efficiency: Time vs Memory\n(Red = Eventalign-based, Orange = F5C Index-based, Blue = Alignment-only)',
                     fontweight='bold')
        ax.grid(alpha=0.3)

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 4: Summary table
        fig, ax = plt.subplots(figsize=(14, 10))
        ax.axis('off')

        # Create summary table
        summary = agg_df.groupby('tool').agg({
            's_total': ['mean', 'std'],
            'cpu_time_total': ['mean', 'std'],
            'max_rss_total': ['mean', 'max'],
            'io_in_total': 'sum',
            'io_out_total': 'sum'
        }).round(2)

        summary.columns = ['_'.join(col).strip() for col in summary.columns.values]
        summary = summary.reset_index()

        # Rename columns for display
        col_names = {
            'tool': 'Tool',
            's_total_mean': 'Time Mean (s)',
            's_total_std': 'Time Std (s)',
            'cpu_time_total_mean': 'CPU Mean (s)',
            'cpu_time_total_std': 'CPU Std (s)',
            'max_rss_total_mean': 'Mem Mean (MB)',
            'max_rss_total_max': 'Mem Max (MB)',
            'io_in_total_sum': 'I/O Read (MB)',
            'io_out_total_sum': 'I/O Write (MB)'
        }
        summary = summary.rename(columns=col_names)

        table = ax.table(cellText=summary.values,
                        colLabels=summary.columns,
                        cellLoc='center',
                        loc='center')

        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.2, 1.5)

        # Style header
        for i in range(len(summary.columns)):
            table[(0, i)].set_facecolor('#3498db')
            table[(0, i)].set_text_props(color='white', fontweight='bold')

        ax.set_title('Resource Summary by Tool (including prerequisites)',
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
        parser = argparse.ArgumentParser(description='Aggregate resources by tool with dependencies')
        parser.add_argument('--benchmark-dir', required=True, help='Directory with benchmark files')
        parser.add_argument('--output-tsv', required=True, help='Output TSV file')
        parser.add_argument('--output-pdf', required=True, help='Output PDF visualization')
        args = parser.parse_args()

        benchmark_dir = args.benchmark_dir
        output_tsv = args.output_tsv
        output_pdf = args.output_pdf

    print(f"Loading benchmark files from {benchmark_dir}")
    df = load_benchmark_files(benchmark_dir)

    if df.empty:
        print("No benchmark data found")
        pd.DataFrame().to_csv(output_tsv, sep='\t', index=False)
        return

    print(f"Loaded {len(df)} benchmark records")
    print(f"Rules found: {df['rule'].unique().tolist()}")

    # Aggregate by tool with dependencies
    print("Aggregating by tool with dependencies...")
    agg_df = aggregate_by_tool_with_dependencies(df)

    if agg_df.empty:
        print("No modification tool data found")
        pd.DataFrame().to_csv(output_tsv, sep='\t', index=False)
        return

    print(f"Aggregated {len(agg_df)} records for {agg_df['tool'].nunique()} tools")

    # Save TSV
    agg_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Saved aggregated data to {output_tsv}")

    # Create visualization
    print("Creating visualization...")
    create_resource_visualization(agg_df, output_pdf)


if __name__ == '__main__':
    main()
