"""
Generate comprehensive PDF report for RNA modification detection benchmarking.

Features:
- Multi-page PDF with overall, per-modification type, and per-tool analysis
- Bar plots comparing tools by F1, precision, recall, AUPRC, AUROC
- Tables with detailed metrics
- Support for multi-window evaluation
- Tool ranking tables
- Optimal threshold analysis

Output:
- benchmark_report.pdf - comprehensive PDF report
"""

import os
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
except ImportError:
    print("Error: matplotlib/seaborn not available. Install with: pip install matplotlib seaborn")
    raise

# PDF page settings
PAGE_WIDTH = 11.69  # A4 landscape width in inches
PAGE_HEIGHT = 8.27  # A4 landscape height in inches


def safe_float(val):
    """Convert value to float, return NaN if invalid."""
    try:
        return float(val) if not pd.isna(val) else np.nan
    except (TypeError, ValueError):
        return np.nan


def load_benchmark_data(benchmark_dir):
    """Load all benchmark data files from directory."""
    data = {}

    # Load accuracy summary (per modification type)
    accuracy_path = os.path.join(benchmark_dir, "accuracy_summary.tsv")
    if os.path.exists(accuracy_path):
        data['accuracy'] = pd.read_csv(accuracy_path, sep='\t')
    else:
        data['accuracy'] = pd.DataFrame()

    # Load overall accuracy
    overall_path = os.path.join(benchmark_dir, "accuracy_summary_overall.tsv")
    if os.path.exists(overall_path):
        data['overall'] = pd.read_csv(overall_path, sep='\t')
    else:
        data['overall'] = pd.DataFrame()

    # Load threshold evaluation
    threshold_path = os.path.join(benchmark_dir, "threshold_evaluation.tsv")
    if os.path.exists(threshold_path):
        data['threshold'] = pd.read_csv(threshold_path, sep='\t')
    else:
        data['threshold'] = pd.DataFrame()

    # Load optimal thresholds
    optimal_path = os.path.join(benchmark_dir, "optimal_thresholds.tsv")
    if os.path.exists(optimal_path):
        data['optimal'] = pd.read_csv(optimal_path, sep='\t')
    else:
        data['optimal'] = pd.DataFrame()

    # Load optimal thresholds detailed
    optimal_detail_path = os.path.join(benchmark_dir, "optimal_thresholds_detailed.tsv")
    if os.path.exists(optimal_detail_path):
        data['optimal_detail'] = pd.read_csv(optimal_detail_path, sep='\t')
    else:
        data['optimal_detail'] = data['optimal']  # Fallback

    # Load score distributions
    dist_path = os.path.join(benchmark_dir, "score_distributions.tsv")
    if os.path.exists(dist_path):
        data['distributions'] = pd.read_csv(dist_path, sep='\t')
    else:
        data['distributions'] = pd.DataFrame()

    # Load optimal score per tool
    optimal_score_path = os.path.join(benchmark_dir, "optimal_score_per_tool.tsv")
    if os.path.exists(optimal_score_path):
        data['optimal_score'] = pd.read_csv(optimal_score_path, sep='\t')
    else:
        data['optimal_score'] = pd.DataFrame()

    # Load all scores evaluation (for score column discovery)
    all_scores_path = os.path.join(benchmark_dir, "all_scores_evaluation.tsv")
    if os.path.exists(all_scores_path):
        data['all_scores'] = pd.read_csv(all_scores_path, sep='\t')
    else:
        data['all_scores'] = pd.DataFrame()

    # Load per-comparison accuracy (for scatterplot by comparison)
    by_comparison_path = os.path.join(benchmark_dir, "accuracy_summary_by_comparison.tsv")
    if os.path.exists(by_comparison_path):
        data['by_comparison'] = pd.read_csv(by_comparison_path, sep='\t')
    else:
        data['by_comparison'] = pd.DataFrame()

    return data


def plot_roc_prc_curves(data, pdf, window=0):
    """Plot ROC and PR curves for all tools at a given window.

    This visualizes the trade-off between precision and recall (PR curve)
    and true positive rate vs false positive rate (ROC curve) for each tool.
    """
    df = data.get('roc_pr_data', pd.DataFrame())

    if df.empty:
        # Create placeholder page explaining data requirements
        fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
        ax = fig.add_subplot(111)
        ax.axis('off')

        ax.text(0.5, 0.5, 'ROC/PR Curves\n\nRequires detailed predictions with scores.\n'
                'Run benchmark_detailed.py to generate required data.',
                fontsize=14, ha='center', va='center', transform=ax.transAxes)

        ax.text(0.5, 0.95, 'ROC and Precision-Recall Curves',
                fontsize=18, fontweight='bold', ha='center', transform=ax.transAxes)

        # Add note about data location
        ax.text(0.5, 0.15, 'Data source: {project}/results/benchmarks/detailed_predictions.tsv',
                fontsize=10, ha='center', va='center', transform=ax.transAxes,
                color='gray', style='italic')

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        return

    # Filter by window if available
    if 'window' in df.columns:
        df = df[df['window'] == window].copy()

    if df.empty:
        return

    tools = df['tool'].unique()

    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(PAGE_WIDTH, PAGE_HEIGHT * 0.8))

    colors = plt.cm.tab10(np.linspace(0, 1, len(tools)))

    # Plot ROC curves
    ax_roc = axes[0]
    for i, tool in enumerate(tools):
        tool_df = df[df['tool'] == tool]

        if 'fpr' in tool_df.columns and 'tpr' in tool_df.columns:
            ax_roc.plot(tool_df['fpr'], tool_df['tpr'],
                       label=f"{tool} (AUC={tool_df['auroc'].iloc[0]:.3f})" if 'auroc' in tool_df.columns else tool,
                       color=colors[i], linewidth=2)

    ax_roc.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random')
    ax_roc.set_xlabel('False Positive Rate', fontsize=11)
    ax_roc.set_ylabel('True Positive Rate', fontsize=11)
    ax_roc.set_title(f'ROC Curves (window={window}nt)', fontsize=14, fontweight='bold')
    ax_roc.legend(loc='lower right', fontsize=9)
    ax_roc.set_xlim(0, 1)
    ax_roc.set_ylim(0, 1)
    ax_roc.grid(alpha=0.3)

    # Plot PR curves
    ax_pr = axes[1]
    for i, tool in enumerate(tools):
        tool_df = df[df['tool'] == tool]

        if 'recall' in tool_df.columns and 'precision' in tool_df.columns:
            ax_pr.plot(tool_df['recall'], tool_df['precision'],
                      label=f"{tool} (AUC={tool_df['auprc'].iloc[0]:.3f})" if 'auprc' in tool_df.columns else tool,
                      color=colors[i], linewidth=2)

    ax_pr.set_xlabel('Recall (Sensitivity)', fontsize=11)
    ax_pr.set_ylabel('Precision', fontsize=11)
    ax_pr.set_title(f'Precision-Recall Curves (window={window}nt)', fontsize=14, fontweight='bold')
    ax_pr.legend(loc='lower left', fontsize=9)
    ax_pr.set_xlim(0, 1)
    ax_pr.set_ylim(0, 1)
    ax_pr.grid(alpha=0.3)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def plot_metrics_by_modification_type(data, pdf, window=0):
    """Plot metrics grouped by modification type as grouped bar chart.

    Data source: {project}/results/benchmarks/accuracy_summary.tsv
    """
    df = data['accuracy']

    if df.empty:
        return

    # Filter by window if available
    if 'window' in df.columns:
        df = df[df['window'] == window].copy()

    if df.empty or 'modification_type' not in df.columns:
        return

    mod_types = df['modification_type'].unique()
    tools = df['tool'].unique()

    if len(mod_types) == 0 or len(tools) == 0:
        return

    # Create figure with multiple subplots for key metrics
    metrics = ['f1', 'precision', 'recall', 'auprc', 'auroc']
    metrics = [m for m in metrics if m in df.columns]

    if not metrics:
        return

    n_metrics = len(metrics)
    fig, axes = plt.subplots(1, n_metrics, figsize=(PAGE_WIDTH, PAGE_HEIGHT * 0.6))

    if n_metrics == 1:
        axes = [axes]

    colors = plt.cm.Set2(np.linspace(0, 1, len(tools)))

    for ax, metric in zip(axes, metrics):
        x = np.arange(len(mod_types))
        width = 0.8 / len(tools)

        for i, tool in enumerate(tools):
            tool_df = df[df['tool'] == tool]
            vals = [tool_df[tool_df['modification_type'] == mt][metric].mean()
                   for mt in mod_types]
            ax.bar(x + i * width, vals, width, label=tool, color=colors[i])

        ax.set_xlabel('Modification Type', fontsize=10)
        ax.set_ylabel(metric.upper(), fontsize=10)
        ax.set_title(f'{metric.upper()} by Modification Type', fontsize=12, fontweight='bold')
        ax.set_xticks(x + width * (len(tools) - 1) / 2)
        ax.set_xticklabels(mod_types, rotation=45, ha='right', fontsize=9)
        ax.set_ylim(0, 1)
        ax.grid(axis='y', alpha=0.3)

    # Add single legend for all subplots
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=min(len(tools), 5),
               bbox_to_anchor=(0.5, 1.02), fontsize=9)

    # Add data source annotation
    fig.text(0.5, -0.02, 'Data source: {project}/results/benchmarks/accuracy_summary.tsv',
             ha='center', fontsize=9, color='gray', style='italic')

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def plot_tool_called_sites_comparison(data, pdf, tool, window=0):
    """Plot called sites for each comparison (e.g., control_0_case_0) for a tool.

    Uses accuracy_summary.tsv data grouped by modification_type as a proxy for comparisons.
    """
    df = data['accuracy']

    if df.empty:
        return

    # Filter for this tool and window
    tool_df = df[df['tool'] == tool].copy()
    if 'window' in tool_df.columns:
        tool_df = tool_df[tool_df['window'] == window]

    if tool_df.empty:
        return

    # Check if called_sites column exists
    if 'called_sites' not in tool_df.columns:
        fig = plt.figure(figsize=(PAGE_WIDTH * 0.9, PAGE_HEIGHT * 0.4))
        ax = fig.add_subplot(111)
        ax.axis('off')
        ax.text(0.5, 0.5, f'{tool}: Called sites data not available',
                ha='center', va='center', fontsize=12, color='gray')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        return

    # Use modification_type as comparison groups
    if 'modification_type' in tool_df.columns:
        groups = tool_df.groupby('modification_type')['called_sites'].sum().sort_values(ascending=True)
    else:
        # Single value
        groups = pd.Series({'overall': tool_df['called_sites'].sum()})

    fig, ax = plt.subplots(figsize=(PAGE_WIDTH * 0.9, PAGE_HEIGHT * 0.4))

    bars = ax.barh(range(len(groups)), groups.values, color='steelblue', edgecolor='navy')

    ax.set_xlabel('Called Sites', fontsize=11)
    ax.set_ylabel('Modification Type', fontsize=11)
    ax.set_title(f'{tool}: Called Sites per Comparison (window={window}nt)', fontsize=14, fontweight='bold')
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups.index, fontsize=9)
    ax.grid(axis='x', alpha=0.3)

    # Add count labels on bars
    for bar, count in zip(bars, groups.values):
        ax.text(bar.get_width() + max(groups.values) * 0.01,
               bar.get_y() + bar.get_height()/2,
               f'{int(count):,}', ha='left', va='center', fontsize=8)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def plot_tool_auprc_auroc_scatter(data, pdf, tool, window=0):
    """Plot AUPRC vs AUROC scatter for each comparison (not modification_type).

    Uses per-comparison metrics (accuracy_summary_by_comparison.tsv) which shows
    how each tool performs on different sample pairs (native_0_control_0, etc).
    """
    # Use per-comparison data for scatterplot
    df = data.get('by_comparison', pd.DataFrame())

    if df.empty:
        # Fallback to accuracy data if by_comparison not available
        df = data.get('accuracy', pd.DataFrame())

    if df.empty:
        return

    tool_df = df[(df['tool'] == tool) & (df['window'] == window)].copy()

    if tool_df.empty:
        return

    if 'auprc' not in tool_df.columns or 'auroc' not in tool_df.columns:
        return

    # Remove NaN values
    tool_df = tool_df.dropna(subset=['auprc', 'auroc'])

    if tool_df.empty:
        return

    fig, ax = plt.subplots(figsize=(PAGE_WIDTH * 0.7, PAGE_HEIGHT * 0.5))

    # Plot scatter with colors by comparison or modification_type
    label_col = None
    color_col = None
    for col in ['comparison', 'modification_type']:
        if col in tool_df.columns:
            label_col = col
            color_col = col
            break

    if color_col and tool_df[color_col].nunique() > 1:
        # Color by comparison/modification_type
        unique_vals = tool_df[color_col].unique()
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_vals)))
        color_map = dict(zip(unique_vals, colors))

        for val in unique_vals:
            subset = tool_df[tool_df[color_col] == val]
            ax.scatter(subset['auroc'], subset['auprc'], s=100, alpha=0.7,
                      c=[color_map[val]], label=str(val)[:20], edgecolors='navy')
        ax.legend(loc='lower right', fontsize=8, title=color_col.replace('_', ' ').title())
    else:
        # Single color
        ax.scatter(tool_df['auroc'], tool_df['auprc'], s=100, alpha=0.7,
                  c='steelblue', edgecolors='navy')

    # Add diagonal reference line
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='AUPRC = AUROC')

    # Add labels for each point
    if label_col:
        for _, row in tool_df.iterrows():
            ax.annotate(str(row[label_col])[:15],
                       (row['auroc'], row['auprc']),
                       xytext=(3, 3), textcoords='offset points',
                       fontsize=8, alpha=0.7)

    ax.set_xlabel('AUROC', fontsize=11)
    ax.set_ylabel('AUPRC', fontsize=11)
    ax.set_title(f'{tool}: AUPRC vs AUROC (by {label_col.replace("_", " ") if label_col else "sample"})',
                fontsize=14, fontweight='bold')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)

    # Add mean point
    mean_auroc = tool_df['auroc'].mean()
    mean_auprc = tool_df['auprc'].mean()
    ax.scatter([mean_auroc], [mean_auprc], s=200, c='red', marker='*',
              edgecolors='darkred', zorder=5, label=f'Mean ({mean_auroc:.3f}, {mean_auprc:.3f})')

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def display_available_score_columns(data, pdf, tool):
    """Display all available score columns that could be used for ROC calculation.

    Uses optimal_score_per_tool.tsv if available, otherwise shows common patterns.
    """
    df = data.get('optimal_score', pd.DataFrame())

    fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT * 0.5))
    ax = fig.add_subplot(111)
    ax.axis('off')

    # Check if we have score column data for this tool
    if not df.empty and 'tool' in df.columns:
        tool_df = df[df['tool'] == tool]

        if not tool_df.empty and 'original_score_column' in tool_df.columns:
            # Show actual score columns found
            score_cols = tool_df['original_score_column'].unique().tolist()

            text = f"Tool: {tool}\n\n"
            text += "Score columns available for ROC calculation:\n\n"
            for col in score_cols:
                text += f"  • {col}\n"

            if 'score_column' in tool_df.columns:
                used = tool_df['score_column'].iloc[0]
                text += f"\nColumn used for this report: {used}"

            ax.text(0.1, 0.5, text, fontsize=11, va='center', transform=ax.transAxes,
                   family='monospace',
                   bbox=dict(boxstyle='round', facecolor='#d4edda', alpha=0.5))
        else:
            # Show common patterns
            text = f"Tool: {tool}\n\n"
            text += "Common score columns for this tool type:\n\n"
            text += "  • p_value / pvalue / pval\n"
            text += "  • adjusted_p_value / padj / FDR\n"
            text += "  • mod_ratio / stoichiometry\n"
            text += "  • diff_mod_rate / mod_diff\n"
            text += "  • z_score / statistic\n"
            text += "  • probability / score\n\n"
            text += "Run benchmark_score_optimization rule\n"
            text += "for detailed score column analysis."

            ax.text(0.1, 0.5, text, fontsize=10, va='center', transform=ax.transAxes,
                   family='monospace',
                   bbox=dict(boxstyle='round', facecolor='#fff3cd', alpha=0.5))
    else:
        # No data available - show generic message
        text = f"Tool: {tool}\n\n"
        text += "Score column analysis requires:\n\n"
        text += "  • Run benchmark_score_optimization rule\n"
        text += "  • Output: optimal_score_per_tool.tsv\n\n"
        text += "Common score columns include:\n"
        text += "  p_value, adjusted_p_value, mod_ratio,\n"
        text += "  stoichiometry, mod_diff, z_score"

        ax.text(0.1, 0.5, text, fontsize=10, va='center', transform=ax.transAxes,
               family='monospace',
               bbox=dict(boxstyle='round', facecolor='#ecf0f1', alpha=0.5))

    ax.set_title(f'{tool}: Available Score Columns for ROC Calculation',
                fontsize=14, fontweight='bold')

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def add_title_page(pdf, data):
    """Add title page to PDF."""
    fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
    ax = fig.add_subplot(111)
    ax.axis('off')

    # Title
    ax.text(0.5, 0.75, 'NanoRNAMod Benchmark Report',
            fontsize=32, fontweight='bold', ha='center', va='center',
            transform=ax.transAxes)

    # Subtitle
    ax.text(0.5, 0.60, 'RNA Modification Detection Tool Comparison',
            fontsize=18, ha='center', va='center', transform=ax.transAxes,
            color='#34495e')

    # Date
    ax.text(0.5, 0.48, f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}',
            fontsize=12, ha='center', va='center', transform=ax.transAxes,
            color='gray')

    # Summary stats box
    summary_lines = []

    if not data['accuracy'].empty:
        tools = data['accuracy']['tool'].nunique()
        summary_lines.append(f"Tools Evaluated: {tools}")

        if 'modification_type' in data['accuracy'].columns:
            mod_types = data['accuracy']['modification_type'].unique()
            summary_lines.append(f"Modification Types: {', '.join(map(str, mod_types))}")

        if 'window' in data['accuracy'].columns:
            windows = sorted(data['accuracy']['window'].unique())
            summary_lines.append(f"Position Windows: {windows}")

    if not data['threshold'].empty:
        n_thresholds = data['threshold']['threshold'].nunique()
        summary_lines.append(f"Thresholds Evaluated: {n_thresholds}")

    if summary_lines:
        summary_text = "\n".join(summary_lines)
        ax.text(0.5, 0.25, summary_text,
                fontsize=11, ha='center', va='center', transform=ax.transAxes,
                family='monospace',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#ecf0f1', alpha=0.8))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


# Tool prerequisites information
TOOL_PREREQUISITES = {
    'xpore': {
        'inputs': [
            'Native BAM (minimap2 alignment)',
            'Control BAM (minimap2 alignment)',
            'Reference FASTA (transcriptome)',
            'Event alignment TSV (f5c eventalign --collapse)'
        ],
        'conda': 'xpore',
        'container': 'nanocompore/xpore',
        'description': 'Differential RNA modification detection using nanopore direct RNA sequencing',
        'output': 'xpore_results.tsv',
        'score_columns': ['pval_CASE_vs_CONTROL', 'z_score_CASE_vs_CONTROL', 'diff_mod_rate_CASE_vs_CONTROL']
    },
    'nanocompore': {
        'inputs': [
            'Native BAM (minimap2 alignment)',
            'Control BAM (minimap2 alignment)',
            'Reference FASTA (transcriptome)',
            'Event alignment TSV (f5c eventalign --collapse)'
        ],
        'conda': 'nanocompore',
        'container': 'nanocompore/nanocompore',
        'description': 'Nanopore RNA modification calling using signal intensity and dwell time',
        'output': 'nanocompore_results.tsv',
        'score_columns': ['GMM_logit_pvalue', 'KS_dwell_pvalue', 'KS_intensity_pvalue', 'Logit_LOR']
    },
    'baleen': {
        'inputs': [
            'Native BAM (minimap2 alignment)',
            'Control BAM (minimap2 alignment)',
            'Reference FASTA (transcriptome)',
            'Event alignment TSV (f5c eventalign --collapse)'
        ],
        'conda': 'baleen',
        'container': ' nejlab/baleen',
        'description': 'RNA modification detection using k-mer based models',
        'output': 'site_results.tsv',
        'score_columns': ['p_value', 'padj', 'effect_size', 'stoichiometry']
    },
    'pybaleen': {
        'inputs': [
            'Native BAM (minimap2 alignment)',
            'Control BAM (minimap2 alignment)',
            'Native FASTQ (basecalled reads)',
            'Control FASTQ (basecalled reads)',
            'Native BLOW5 (raw signal)',
            'Control BLOW5 (raw signal)',
            'Reference FASTA (transcriptome)'
        ],
        'conda': 'pybaleen',
        'container': 'pybaleen',
        'description': 'CUDA-accelerated DTW + HMM modification detection',
        'output': 'site_results.tsv',
        'score_columns': ['mod_ratio', 'pvalue', 'padj', 'effect_size', 'mean_p_mod', 'stoichiometry']
    },
    'differr': {
        'inputs': [
            'Native BAM (minimap2 alignment)',
            'Control BAM (minimap2 alignment)',
            'Reference FASTA (transcriptome)'
        ],
        'conda': 'differr',
        'container': '医用',
        'description': 'Differential error rate analysis for RNA modification detection',
        'output': 'differr_results.tsv',
        'score_columns': ['-log10 P value', '-log10 FDR', 'odds ratio', 'G statistic']
    },
    'drummer': {
        'inputs': [
            'Native BAM (minimap2 alignment)',
            'Control BAM (minimap2 alignment)',
            'Reference FASTA (transcriptome)'
        ],
        'conda': 'drummer',
        'container': '暂无',
        'description': 'DRUMMER - Differential RNA Modification detection',
        'output': 'drummer_results.tsv',
        'score_columns': ['OR_padj', 'G_padj', 'G_test', 'log2_(OR)']
    },
    'eligos2': {
        'inputs': [
            'Native BAM (minimap2 alignment)',
            'Control BAM (minimap2 alignment)',
            'Reference FASTA (transcriptome)',
            'Feature annotations (GTF)'
        ],
        'conda': 'eligos2',
        'container': 'nanopore/eligos2',
        'description': 'Epitranscriptional landscape inferring from RNA sequencing',
        'output': 'eligos2_results.tsv',
        'score_columns': ['p_value', 'adjusted_p_value', 'mod_ratio', 'coverage']
    },
    'epinano': {
        'inputs': [
            'Native BAM (minimap2 alignment)',
            'Control BAM (minimap2 alignment)',
            'Reference FASTA (transcriptome)'
        ],
        'conda': 'epinano',
        'container': ' Nanopore/EpiNano',
        'description': 'Detection of m6A RNA modifications using nanopore direct RNA sequencing',
        'output': 'epinano_results.tsv',
        'score_columns': ['SVM_score', 'probability', 'k_mer_score']
    },
    'tandemmod': {
        'inputs': [
            'Native BAM (minimap2 alignment)',
            'Native FASTQ (basecalled reads)',
            'Reference FASTA (transcriptome)'
        ],
        'conda': 'tandemmod',
        'container': ' psmsl/tandemmod',
        'description': 'RNA modification detection without control sample',
        'output': 'tandemmod_results.tsv',
        'score_columns': ['mod_score', 'mod_prob', 'coverage']
    },
    'directrm': {
        'inputs': [
            'Native BAM (minimap2 alignment)',
            'Native FASTQ (basecalled reads)',
            'Reference FASTA (transcriptome)'
        ],
        'conda': 'directrm',
        'container': '暂无',
        'description': 'Direct RNA modification detection',
        'output': 'directrm_results.tsv',
        'score_columns': ['score', 'probability', 'coverage']
    },
    'm6atm': {
        'inputs': [
            'Native BAM (minimap2 alignment)',
            'Native FASTQ (basecalled reads)',
            'Reference FASTA (transcriptome)'
        ],
        'conda': 'm6atm',
        'container': '暂无',
        'description': 'm6A detection using transformer model',
        'output': 'm6atm_results.tsv',
        'score_columns': ['score', 'probability']
    },
    'rnano': {
        'inputs': [
            'Native BAM (minimap2 alignment)',
            'Native FASTQ (basecalled reads)',
            'Reference FASTA (transcriptome)'
        ],
        'conda': 'rnano',
        'container': '暂无',
        'description': 'RNA modification detection using RNN',
        'output': 'rnano_results.tsv',
        'score_columns': ['score', 'probability', 'coverage']
    }
}


def add_tool_prerequisites_page(pdf, data):
    """Add page showing prerequisites for each tool."""
    fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
    ax = fig.add_subplot(111)
    ax.axis('off')

    ax.text(0.5, 0.96, 'Tool Prerequisites',
            fontsize=20, fontweight='bold', ha='center', va='center',
            transform=ax.transAxes)

    ax.text(0.5, 0.92, 'Input requirements for each modification detection tool',
            fontsize=11, ha='center', va='center', transform=ax.transAxes,
            color='gray')

    # Get tools from data
    df = data['accuracy']
    if df.empty:
        ax.text(0.5, 0.5, 'No tool data available',
                fontsize=12, ha='center', va='center', transform=ax.transAxes)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        return

    tools = sorted(df['tool'].unique())

    # Create table data
    table_data = []
    headers = ['Tool', 'Input Files', 'Score Columns', 'Output']

    for tool in tools:
        info = TOOL_PREREQUISITES.get(tool, {})
        inputs = info.get('inputs', ['Unknown'])
        score_cols = info.get('score_columns', ['Unknown'])
        output = info.get('output', 'Unknown')

        # Truncate for display
        inputs_str = '\n'.join(inputs[:3])
        if len(inputs) > 3:
            inputs_str += f'\n+{len(inputs)-3} more'
        scores_str = ', '.join(score_cols[:3])
        if len(score_cols) > 3:
            scores_str += f' +{len(score_cols)-3}'

        table_data.append([tool, inputs_str, scores_str, output])

    # Create table
    table = ax.table(cellText=table_data,
                     colLabels=headers,
                     cellLoc='left',
                     loc='center',
                     colWidths=[0.12, 0.45, 0.28, 0.15])

    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.1, 2.0)

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

    ax.set_title('', fontsize=14, fontweight='bold')

    # Add legend for common input types
    legend_text = """Common Input Types:
• BAM: Aligned reads (minimap2)
• FASTQ: Basecalled reads
• BLOW5: Raw nanopore signal
• Event TSV: f5c eventalign --collapse output
• Reference FASTA: Transcriptome sequence
• GTF: Gene annotations"""

    ax.text(0.02, 0.08, legend_text,
            fontsize=8, va='bottom', ha='left',
            transform=ax.transAxes, family='monospace',
            bbox=dict(boxstyle='round', facecolor='#f8f9fa', alpha=0.8))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def add_table_of_contents(pdf):
    """Add table of contents page."""
    fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
    ax = fig.add_subplot(111)
    ax.axis('off')

    ax.text(0.5, 0.90, 'Table of Contents',
            fontsize=20, fontweight='bold', ha='center', va='center',
            transform=ax.transAxes)

    toc_items = [
        "1. Metrics Explanation",
        "2. Tool Prerequisites",
        "3. Overall Summary",
        "4. Tool Comparison Charts",
        "5. Called Sites Comparison",
        "6. ROC/PR Curves Comparison",
        "7. Metrics by Window",
        "8. F1 Score Heatmap",
        "9. AUPRC Heatmap",
        "10. AUROC Heatmap",
        "11. Metrics by Modification Type",
        "12. Tool Rankings",
        "13. Per-Tool Detailed Analysis",
        "14. Optimal Thresholds & Score Columns",
        "15. Score Distributions",
        "16. Data Sources"
    ]

    y_pos = 0.75
    for item in toc_items:
        ax.text(0.3, y_pos, item, fontsize=12, transform=ax.transAxes)
        y_pos -= 0.07

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def add_metrics_explanation_page(pdf):
    """Add page explaining all metrics used in the report."""
    fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
    ax = fig.add_subplot(111)
    ax.axis('off')

    ax.text(0.5, 0.95, 'Metrics Explanation',
            fontsize=20, fontweight='bold', ha='center', va='center',
            transform=ax.transAxes)

    explanations = """
Classification Metrics (range 0-1, higher is better):

  • Precision (Positive Predictive Value)
    TP / (TP + FP) - Fraction of predicted sites that are true modifications.
    High precision = few false positives.

  • Recall (Sensitivity, True Positive Rate)
    TP / (TP + FN) - Fraction of true modifications that were detected.
    High recall = few false negatives.

  • F1 Score
    2 × (Precision × Recall) / (Precision + Recall) - Harmonic mean of precision
    and recall. Balances both metrics; useful when classes are imbalanced.

  • AUPRC (Area Under Precision-Recall Curve)
    Integral of precision vs recall across all thresholds. Better than AUROC
    for imbalanced datasets where positive sites are rare.

  • AUROC (Area Under ROC Curve)
    Area under TPR vs FPR curve. Probability that a random positive site
    scores higher than a random negative site.

  • MCC (Matthews Correlation Coefficient)
    Correlation coefficient between predicted and actual classifications.
    Range -1 to +1. More robust for imbalanced data than accuracy.
    NOTE: MCC is shown as "not available" when ground truth does not contain
    explicit negative sites. In this case, TN cannot be computed accurately.


Site Counts:

  • Called Sites
    Total number of modification sites reported by a tool. Important for
    comparing tools - a tool calling 1000 sites vs one calling 10 sites
    may have different precision/recall trade-offs even with similar scores.

  • Total Truth Sites
    Number of known modification sites in the ground truth.


Positional Tolerance:

  • Window (nt)
    Positional tolerance in nucleotides. A predicted site at position P matches
    a true site at position T if |P - T| ≤ window. window=0 requires exact match.


Tool-Specific Scores:

  • Score Column
    Each tool outputs its own scoring metric (p-value, probability, effect size,
    etc.). The benchmark uses these native scores for thresholding. Lower p-values
    or higher probabilities indicate stronger modification signals.

  • Threshold
    The score cutoff used to call a site as modified. Optimal thresholds are
    determined by maximizing F1 score across all tested thresholds.
"""

    ax.text(0.05, 0.85, explanations,
            fontsize=9, family='monospace', va='top',
            transform=ax.transAxes,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8f9fa', alpha=0.8))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def plot_overall_metrics_bar(data, pdf, window=0):
    """Create bar plot comparing all tools on overall metrics."""
    # Use overall data if available, otherwise aggregate accuracy
    df = data['overall'] if not data['overall'].empty else data['accuracy']

    if df.empty:
        return

    # Filter by window if available
    if 'window' in df.columns:
        df = df[df['window'] == window].copy()

    if df.empty:
        return

    # Aggregate by tool — only include columns that actually exist
    agg_spec = {col: 'mean' for col in ['precision', 'recall', 'f1', 'auprc', 'auroc', 'mcc']
                if col in df.columns}
    if not agg_spec:
        return
    agg_df = df.groupby('tool').agg(agg_spec).reset_index()

    fig, axes = plt.subplots(2, 3, figsize=(PAGE_WIDTH, PAGE_HEIGHT))
    fig.suptitle(f'Overall Tool Comparison (window={window}nt)', fontsize=16, fontweight='bold')

    metrics = [
        ('f1', 'F1 Score', axes[0, 0]),
        ('precision', 'Precision', axes[0, 1]),
        ('recall', 'Recall', axes[0, 2]),
        ('auprc', 'AUPRC', axes[1, 0]),
        ('auroc', 'AUROC', axes[1, 1]),
        ('mcc', 'MCC', axes[1, 2])
    ]

    tools = sorted(agg_df['tool'].unique())
    colors = plt.cm.Set2(np.linspace(0, 1, max(len(tools), 1)))
    color_map = dict(zip(tools, colors))

    for metric, title, ax in metrics:
        if metric not in agg_df.columns or agg_df[metric].isna().all():
            ax.axis('off')
            ax.text(0.5, 0.5, f'{title}\n(not available)', ha='center', va='center',
                   fontsize=10, color='gray')
            continue

        plot_df = agg_df.dropna(subset=[metric]).sort_values(metric, ascending=True)

        if plot_df.empty:
            ax.axis('off')
            continue

        bars = ax.barh(plot_df['tool'], plot_df[metric],
                       color=[color_map.get(t, 'gray') for t in plot_df['tool']])
        ax.set_xlabel(title, fontsize=10)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xlim(0, 1.0)
        ax.grid(axis='x', alpha=0.3)

        # Add value labels
        for bar, val in zip(bars, plot_df[metric]):
            if not np.isnan(val):
                ax.text(val + 0.02, bar.get_y() + bar.get_height()/2,
                       f'{val:.3f}', va='center', fontsize=8)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def plot_called_sites_comparison(data, pdf, window=0):
    """Create bar plot comparing total called sites per tool."""
    df = data['overall'] if not data['overall'].empty else data['accuracy']

    if df.empty:
        return

    # Filter by window if available
    if 'window' in df.columns:
        df = df[df['window'] == window].copy()

    if df.empty:
        return

    # Check if called_sites column exists
    if 'called_sites' not in df.columns:
        return

    # Aggregate by tool
    agg_df = df.groupby('tool').agg({
        'called_sites': 'sum',
        'total_truth': 'first' if 'total_truth' in df.columns else 'sum'
    }).reset_index()

    fig, axes = plt.subplots(1, 2, figsize=(PAGE_WIDTH, PAGE_HEIGHT * 0.5))
    fig.suptitle(f'Site Counts Comparison (window={window}nt)', fontsize=14, fontweight='bold')

    tools = sorted(agg_df['tool'].unique())
    colors = plt.cm.Set2(np.linspace(0, 1, max(len(tools), 1)))
    color_map = dict(zip(tools, colors))

    # Plot 1: Called sites
    ax = axes[0]
    plot_df = agg_df.sort_values('called_sites', ascending=True)
    bars = ax.barh(plot_df['tool'], plot_df['called_sites'],
                   color=[color_map.get(t, 'gray') for t in plot_df['tool']])
    ax.set_xlabel('Number of Sites', fontsize=10)
    ax.set_title('Total Called Sites per Tool', fontsize=12, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)

    # Add value labels
    for bar, val in zip(bars, plot_df['called_sites']):
        ax.text(val + max(plot_df['called_sites']) * 0.01, bar.get_y() + bar.get_height()/2,
               f'{int(val):,}', va='center', fontsize=8)

    # Plot 2: Called sites vs truth sites
    ax = axes[1]
    if 'total_truth' in agg_df.columns and not agg_df['total_truth'].isna().all():
        truth_sites = agg_df['total_truth'].iloc[0]  # Same for all tools
        ax.axvline(truth_sites, color='red', linestyle='--', linewidth=2, label=f'Truth Sites: {int(truth_sites):,}')
        ax.barh(plot_df['tool'], plot_df['called_sites'],
                color=[color_map.get(t, 'gray') for t in plot_df['tool']], alpha=0.7)
        ax.set_xlabel('Number of Sites', fontsize=10)
        ax.set_title('Called Sites vs Truth Sites', fontsize=12, fontweight='bold')
        ax.legend(loc='lower right', fontsize=8)
        ax.grid(axis='x', alpha=0.3)
    else:
        ax.axis('off')
        ax.text(0.5, 0.5, 'Truth sites count\nnot available', ha='center', va='center',
               fontsize=10, color='gray')

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def add_data_sources_page(pdf, benchmark_dir):
    """Add page listing all data source files."""
    fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
    ax = fig.add_subplot(111)
    ax.axis('off')

    ax.text(0.5, 0.95, 'Data Sources',
            fontsize=20, fontweight='bold', ha='center', va='center',
            transform=ax.transAxes)

    source_text = f"""
All data used to generate this report is stored in the benchmark directory:

  {benchmark_dir}/

Key data files:

  • accuracy_summary.tsv
    Per-modification type metrics for each tool at each window tolerance.

  • accuracy_summary_overall.tsv
    Aggregated overall metrics for each tool at each window tolerance.

  • threshold_evaluation.tsv
    Detailed metrics (P/R/F1/AUPRC/AUROC) computed at multiple score thresholds
    for each tool and modification type.

  • optimal_thresholds.tsv
    Optimal score thresholds determined by maximizing F1 score.

  • optimal_thresholds_detailed.tsv
    Detailed optimal threshold information including original score columns.

  • score_distributions.tsv
    Statistics (min/max/mean/std/median) of score distributions for each tool.

  • resource_by_tool.tsv
    Computational resource usage (time, memory, I/O) aggregated by tool.


Report Generation:

  • Report generated by: benchmark_pdf_report.py
  • Date: {datetime.now().strftime("%Y-%m-%d %H:%M")}


Reproducing the Analysis:

  To reproduce this analysis or access the raw data:

  1. Navigate to the benchmark directory
  2. Load TSV files with pandas: pd.read_csv(file, sep='\\t')
  3. Each file contains tool, modification_type, window columns for filtering
"""

    ax.text(0.05, 0.82, source_text,
            fontsize=9, family='monospace', va='top',
            transform=ax.transAxes,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8f9fa', alpha=0.8))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def plot_metrics_by_window(data, pdf):
    """Create line plot showing metrics across different windows."""
    df = data['overall'] if not data['overall'].empty else data['accuracy']

    if df.empty or 'window' not in df.columns:
        return

    windows = sorted(df['window'].unique())
    if len(windows) <= 1:
        return

    fig, axes = plt.subplots(2, 2, figsize=(PAGE_WIDTH, PAGE_HEIGHT))
    fig.suptitle('Metrics by Position Tolerance Window', fontsize=16, fontweight='bold')

    metrics = [
        ('f1', 'F1 Score', axes[0, 0]),
        ('precision', 'Precision', axes[0, 1]),
        ('recall', 'Recall', axes[1, 0]),
        ('auprc', 'AUPRC', axes[1, 1])
    ]

    tools = sorted(df['tool'].unique())
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h', 'H', '+', 'x', 'X', 'd']

    for metric, title, ax in metrics:
        if metric not in df.columns:
            ax.axis('off')
            continue

        for i, tool in enumerate(tools):
            tool_df = df[df['tool'] == tool].sort_values('window')
            if tool_df.empty:
                continue
            ax.plot(tool_df['window'], tool_df[metric],
                   marker=markers[i % len(markers)],
                   label=tool, linewidth=2, markersize=8, alpha=0.8)

        ax.set_xlabel('Window (nt)', fontsize=10)
        ax.set_ylabel(title, fontsize=10)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_ylim(0, 1.0)
        ax.grid(alpha=0.3)

        # Only show legend on first plot
        if metric == 'f1':
            ax.legend(loc='lower right', fontsize=7, ncol=2)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def plot_heatmap(data, pdf, window=0, metric='f1', vmin=None, vmax=None):
    """Create heatmap of metrics per tool and modification type.

    Parameters:
        vmin, vmax: Override colorbar limits. If None, use data-driven range.
                   Set to 0, 1 for percentage metrics where 0-1 is meaningful.
    """
    df = data['accuracy']

    if df.empty:
        return

    # Filter by window if available
    if 'window' in df.columns:
        df = df[df['window'] == window].copy()

    if df.empty or metric not in df.columns:
        return

    # Create pivot table
    if 'modification_type' in df.columns:
        pivot = df.pivot_table(index='tool', columns='modification_type',
                               values=metric, aggfunc='mean')
        xlabel = 'Modification Type'
    else:
        # Use multiple metrics
        metrics = ['precision', 'recall', 'f1']
        metrics = [m for m in metrics if m in df.columns]
        pivot = df.set_index('tool')[metrics]
        xlabel = 'Metric'

    if pivot.empty:
        return

    # Determine colorbar limits based on data if not specified
    data_min = pivot.min().min()
    data_max = pivot.max().max()

    if vmin is None:
        vmin = data_min
    if vmax is None:
        vmax = data_max

    # For metrics that are inherently 0-1 (precision, recall, f1, etc.),
    # optionally clip to reasonable range
    if metric.lower() in ['precision', 'recall', 'f1', 'auprc', 'auroc', 'specificity']:
        # Use data range but don't extend beyond 0-1
        vmin = max(0, data_min)
        vmax = min(1, data_max)

    fig, ax = plt.subplots(figsize=(PAGE_WIDTH, PAGE_HEIGHT * 0.7))

    # Use reversed colormap: green for low values, red for high values
    # RdYlGn_r = green (low) -> yellow (mid) -> red (high)
    sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn_r',
                vmin=vmin, vmax=vmax, ax=ax,
                cbar_kws={'label': f'{metric.upper()} ({vmin:.2f}-{vmax:.2f})'},
                linewidths=0.5, linecolor='white',
                annot_kws={'fontsize': 10})

    ax.set_title(f'{metric.upper()} by Tool and Modification Type (window={window}nt)',
                 fontsize=14, fontweight='bold')
    ax.set_ylabel('Tool', fontsize=11)
    ax.set_xlabel(xlabel, fontsize=11)

    # Rotate x labels if many columns
    if len(pivot.columns) > 4:
        plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def add_ranking_page(data, pdf, window=0):
    """Add page with tool rankings."""
    df = data['overall'] if not data['overall'].empty else data['accuracy']

    if df.empty:
        return

    # Filter by window
    if 'window' in df.columns:
        df = df[df['window'] == window].copy()

    if df.empty:
        return

    fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
    ax = fig.add_subplot(111)
    ax.axis('off')

    # Compute rankings (1 = best)
    metrics = ['f1', 'precision', 'recall', 'auprc', 'auroc']
    metrics = [m for m in metrics if m in df.columns]

    if not metrics:
        return

    rankings = {}
    for metric in metrics:
        # Higher is better for all metrics, so rank descending
        ranked = df.groupby('tool')[metric].mean().rank(ascending=False)
        rankings[metric] = ranked

    rank_df = pd.DataFrame(rankings)

    # Add average rank
    rank_df['Avg Rank'] = rank_df.mean(axis=1)
    rank_df = rank_df.sort_values('Avg Rank')

    # Rename columns for display
    rank_df.columns = [c.upper() if c != 'Avg Rank' else c for c in rank_df.columns]

    # Reset index to get tool as column
    rank_df = rank_df.reset_index()

    # Create table
    table = ax.table(cellText=rank_df.round(1).values,
                     colLabels=['Tool'] + list(rank_df.columns[1:]),
                     cellLoc='center',
                     loc='center')

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.3, 1.8)

    # Style header
    for i in range(len(rank_df.columns)):
        table[(0, i)].set_facecolor('#3498db')
        table[(0, i)].set_text_props(color='white', fontweight='bold')

    # Highlight top 3
    for i in range(1, min(4, len(rank_df) + 1)):
        for j in range(len(rank_df.columns)):
            if i == 1:
                table[(i, j)].set_facecolor('#d4edda')  # Gold
            elif i == 2:
                table[(i, j)].set_facecolor('#fff3cd')  # Silver
            elif i == 3:
                table[(i, j)].set_facecolor('#cce5ff')  # Bronze

    ax.set_title(f'Tool Rankings (window={window}nt)\nLower rank = better performance',
                 fontsize=16, fontweight='bold', y=0.95)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def add_metrics_table_page(pdf, title, df, window=None, max_rows=30):
    """Add a page with metrics table."""
    if df.empty:
        return

    # Filter by window if specified
    if window is not None and 'window' in df.columns:
        df = df[df['window'] == window].copy()

    if df.empty:
        return

    fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
    ax = fig.add_subplot(111)
    ax.axis('off')

    # Prepare table columns
    display_cols = ['tool']
    if 'modification_type' in df.columns:
        display_cols.append('modification_type')
    display_cols.extend(['precision', 'recall', 'f1', 'auprc', 'auroc'])

    # Add called_sites if available
    if 'called_sites' in df.columns:
        display_cols.append('called_sites')

    # Filter to existing columns
    display_cols = [c for c in display_cols if c in df.columns]

    # Sort by F1 descending
    if 'f1' in df.columns:
        df = df.sort_values('f1', ascending=False)

    table_df = df[display_cols].head(max_rows).copy()

    # Round numeric columns
    for col in ['precision', 'recall', 'f1', 'auprc', 'auroc']:
        if col in table_df.columns:
            table_df[col] = table_df[col].apply(
                lambda x: f'{x:.3f}' if pd.notna(x) else '-'
            )

    # Format called_sites as integer
    if 'called_sites' in table_df.columns:
        table_df['called_sites'] = table_df['called_sites'].apply(
            lambda x: f'{int(x):,}' if pd.notna(x) else '-'
        )

    # Create table
    table = ax.table(cellText=table_df.values,
                     colLabels=display_cols,
                     cellLoc='center',
                     loc='center')

    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.4)

    # Style header
    for i, col in enumerate(display_cols):
        table[(0, i)].set_facecolor('#3498db')
        table[(0, i)].set_text_props(color='white', fontweight='bold')

    # Alternate row colors
    for i in range(1, len(table_df) + 1):
        for j in range(len(display_cols)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#f8f9fa')

    subtitle = f" (window={window}nt)" if window is not None else ""
    ax.set_title(f'{title}{subtitle}', fontsize=14, fontweight='bold', y=0.98)

    if len(df) > max_rows:
        ax.text(0.5, 0.02, f'Showing top {max_rows} of {len(df)} rows',
               fontsize=9, ha='center', transform=ax.transAxes, color='gray')

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def add_per_tool_pages(data, pdf):
    """Add detailed pages for each tool with enhanced visualizations.

    Each tool gets multiple pages:
    - Page 1: Called sites per comparison + AUPRC vs AUROC scatter
    - Page 2: Available score columns for ROC calculation
    """
    df = data['accuracy']

    if df.empty:
        return

    tools = sorted(df['tool'].unique())
    windows = [0]
    if 'window' in df.columns:
        windows = sorted(df['window'].unique())

    for tool in tools:
        tool_df = df[df['tool'] == tool].copy()

        if tool_df.empty:
            continue

        # Page 1: Called sites and AUPRC/AUROC scatter
        fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
        gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

        fig.suptitle(f'Tool: {tool}', fontsize=18, fontweight='bold')

        # Average metrics bar chart (top-left)
        ax_bar = fig.add_subplot(gs[0, 0])

        metrics = ['precision', 'recall', 'f1', 'auprc', 'auroc']
        metrics = [m for m in metrics if m in tool_df.columns]

        if metrics:
            means = tool_df[metrics].mean()
            colors = ['#3498db', '#2ecc71', '#e74c3c', '#9b59b6', '#f39c12']
            bars = ax_bar.bar(metrics, means.values, color=colors[:len(metrics)])
            ax_bar.set_ylim(0, 1.0)
            ax_bar.set_ylabel('Score')
            ax_bar.set_title('Average Metrics', fontsize=12, fontweight='bold')
            ax_bar.tick_params(axis='x', rotation=45)

            for bar, val in zip(bars, means.values):
                if not np.isnan(val):
                    ax_bar.text(bar.get_x() + bar.get_width()/2, val + 0.02,
                               f'{val:.3f}', ha='center', fontsize=9)

        # Called sites by modification type/comparison (top-right)
        ax_sites = fig.add_subplot(gs[0, 1])

        if 'called_sites' in tool_df.columns:
            if 'modification_type' in tool_df.columns:
                # Group by modification type
                sites_by_mod = tool_df.groupby('modification_type')['called_sites'].sum().sort_values(ascending=True)
                bars = ax_sites.barh(sites_by_mod.index, sites_by_mod.values, color='steelblue')
                ax_sites.set_xlabel('Called Sites')
                ax_sites.set_title('Called Sites by Modification Type', fontsize=12, fontweight='bold')

                for bar, val in zip(bars, sites_by_mod.values):
                    ax_sites.text(val + max(sites_by_mod.values) * 0.01, bar.get_y() + bar.get_height()/2,
                                 f'{int(val):,}', va='center', fontsize=8)
            else:
                # Show total
                total = tool_df['called_sites'].sum()
                ax_sites.bar(['Total'], [total], color='steelblue')
                ax_sites.text(0, total + total * 0.02, f'{int(total):,}', ha='center', fontsize=10)
                ax_sites.set_ylabel('Called Sites')
                ax_sites.set_title('Total Called Sites', fontsize=12, fontweight='bold')

            ax_sites.grid(axis='x', alpha=0.3)
        else:
            ax_sites.axis('off')
            ax_sites.text(0.5, 0.5, 'Called sites data\nnot available', ha='center', va='center',
                         fontsize=10, color='gray')

        # AUPRC vs AUROC scatter (bottom-left)
        ax_scatter = fig.add_subplot(gs[1, 0])

        if 'auprc' in tool_df.columns and 'auroc' in tool_df.columns:
            ax_scatter.scatter(tool_df['auroc'], tool_df['auprc'], s=100, alpha=0.7,
                              c='steelblue', edgecolors='navy')
            ax_scatter.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='AUPRC = AUROC')

            # Add labels for each point
            label_col = None
            for col in ['modification_type', 'comparison', 'sample']:
                if col in tool_df.columns:
                    label_col = col
                    break

            if label_col:
                for _, row in tool_df.iterrows():
                    ax_scatter.annotate(str(row[label_col])[:15],
                                       (row['auroc'], row['auprc']),
                                       xytext=(3, 3), textcoords='offset points',
                                       fontsize=7, alpha=0.7)

            # Add mean point
            mean_auroc = tool_df['auroc'].mean()
            mean_auprc = tool_df['auprc'].mean()
            ax_scatter.scatter([mean_auroc], [mean_auprc], s=200, c='red', marker='*',
                              edgecolors='darkred', zorder=5, label=f'Mean ({mean_auroc:.3f}, {mean_auprc:.3f})')

            ax_scatter.set_xlabel('AUROC', fontsize=10)
            ax_scatter.set_ylabel('AUPRC', fontsize=10)
            ax_scatter.set_title('AUPRC vs AUROC', fontsize=12, fontweight='bold')
            ax_scatter.set_xlim(0, 1)
            ax_scatter.set_ylim(0, 1)
            ax_scatter.grid(alpha=0.3)
            ax_scatter.legend(loc='lower right', fontsize=8)
        else:
            ax_scatter.axis('off')
            ax_scatter.text(0.5, 0.5, 'AUPRC/AUROC data\nnot available', ha='center', va='center',
                           fontsize=10, color='gray')

        # Metrics by modification type (bottom-right)
        if 'modification_type' in tool_df.columns:
            ax_mod = fig.add_subplot(gs[1, 1])

            mod_types = tool_df['modification_type'].unique()
            x = np.arange(len(mod_types))
            width = 0.25

            for i, metric in enumerate(['precision', 'recall', 'f1']):
                if metric in tool_df.columns:
                    vals = [tool_df[tool_df['modification_type'] == m][metric].mean()
                           for m in mod_types]
                    ax_mod.bar(x + i * width, vals, width, label=metric)

            ax_mod.set_xlabel('Modification Type')
            ax_mod.set_ylabel('Score')
            ax_mod.set_title('Metrics by Modification Type', fontsize=12, fontweight='bold')
            ax_mod.set_xticks(x + width)
            ax_mod.set_xticklabels(mod_types, rotation=45, ha='right')
            ax_mod.legend(fontsize=8)
            ax_mod.set_ylim(0, 1.0)
            ax_mod.grid(axis='y', alpha=0.3)
        else:
            ax_blank = fig.add_subplot(gs[1, 1])
            ax_blank.axis('off')

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 2: Available score columns
        display_available_score_columns(data, pdf, tool)


def add_optimal_thresholds_page(data, pdf):
    """Add page with optimal thresholds analysis."""
    df = data['optimal_detail'] if not data['optimal_detail'].empty else data['optimal']

    if df.empty:
        return

    fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
    ax = fig.add_subplot(111)
    ax.axis('off')

    # Prepare display columns
    display_cols = ['tool']
    if 'modification_type' in df.columns:
        display_cols.append('modification_type')
    if 'window' in df.columns:
        display_cols.append('window')

    display_cols.extend(['threshold', 'f1', 'precision', 'recall'])

    if 'original_score_column' in df.columns:
        display_cols.append('original_score_column')
    if 'original_threshold' in df.columns:
        display_cols.append('original_threshold')

    # Filter to existing columns
    display_cols = [c for c in display_cols if c in df.columns]

    # Sort by F1
    if 'f1' in df.columns:
        df = df.sort_values('f1', ascending=False)

    table_df = df[display_cols].head(25).copy()

    # Format numeric columns
    for col in ['threshold', 'original_threshold']:
        if col in table_df.columns:
            table_df[col] = table_df[col].apply(
                lambda x: f'{x:.4f}' if pd.notna(x) else '-'
            )
    for col in ['f1', 'precision', 'recall']:
        if col in table_df.columns:
            table_df[col] = table_df[col].apply(
                lambda x: f'{x:.3f}' if pd.notna(x) else '-'
            )

    # Create table
    table = ax.table(cellText=table_df.values,
                     colLabels=display_cols,
                     cellLoc='center',
                     loc='center')

    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1.1, 1.3)

    # Style header
    for i, col in enumerate(display_cols):
        table[(0, i)].set_facecolor('#3498db')
        table[(0, i)].set_text_props(color='white', fontweight='bold')

    ax.set_title(f'Optimal Thresholds by Tool (Top 25 by F1)',
                 fontsize=14, fontweight='bold', y=0.98)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def add_score_distributions_page(data, pdf):
    """Add page with score distribution analysis."""
    df = data['distributions']

    if df.empty:
        return

    fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
    ax = fig.add_subplot(111)
    ax.axis('off')

    # Prepare display columns
    display_cols = ['tool']
    if 'modification_type' in df.columns:
        display_cols.append('modification_type')
    if 'score_column' in df.columns:
        display_cols.append('score_column')
    if 'original_score_column' in df.columns:
        display_cols.append('original_score_column')

    display_cols.extend(['count', 'min', 'max', 'mean', 'std', 'median'])

    # Filter to existing columns
    display_cols = [c for c in display_cols if c in df.columns]

    table_df = df[display_cols].copy()

    # Format numeric columns
    for col in ['min', 'max', 'mean', 'std', 'median']:
        if col in table_df.columns:
            table_df[col] = table_df[col].apply(
                lambda x: f'{x:.4f}' if pd.notna(x) else '-'
            )
    if 'count' in table_df.columns:
        table_df['count'] = table_df['count'].apply(
            lambda x: f'{int(x):,}' if pd.notna(x) else '-'
        )

    # Create table
    table = ax.table(cellText=table_df.values,
                     colLabels=display_cols,
                     cellLoc='center',
                     loc='center')

    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1.1, 1.3)

    # Style header
    for i, col in enumerate(display_cols):
        table[(0, i)].set_facecolor('#3498db')
        table[(0, i)].set_text_props(color='white', fontweight='bold')

    ax.set_title('Score Distributions by Tool',
                 fontsize=14, fontweight='bold', y=0.98)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def generate_pdf_report(data, output_path, benchmark_dir=None):
    """Generate complete PDF report."""
    if benchmark_dir is None:
        benchmark_dir = os.path.dirname(output_path)

    with PdfPages(output_path) as pdf:
        # Title page
        add_title_page(pdf, data)

        # Table of contents
        add_table_of_contents(pdf)

        # Metrics explanation
        add_metrics_explanation_page(pdf)

        # Tool prerequisites (at the beginning for easy checking)
        add_tool_prerequisites_page(pdf, data)

        # Get windows to analyze
        windows = [0]
        if not data['accuracy'].empty and 'window' in data['accuracy'].columns:
            windows = sorted(data['accuracy']['window'].unique())

        # Overall metrics comparison
        for window in windows:
            plot_overall_metrics_bar(data, pdf, window)

        # Called sites comparison
        for window in windows:
            plot_called_sites_comparison(data, pdf, window)

        # ROC/PR Curves comparison
        for window in windows:
            plot_roc_prc_curves(data, pdf, window)

        # Metrics by window (line plot)
        plot_metrics_by_window(data, pdf)

        # Heatmaps - F1, AUPRC, AUROC
        for window in windows:
            plot_heatmap(data, pdf, window, 'f1')

        for window in windows:
            if 'auprc' in data['accuracy'].columns:
                plot_heatmap(data, pdf, window, 'auprc')

        for window in windows:
            if 'auroc' in data['accuracy'].columns:
                plot_heatmap(data, pdf, window, 'auroc')

        # Ranking page
        for window in windows:
            add_ranking_page(data, pdf, window)

        # Overall metrics table
        add_metrics_table_page(pdf, 'Overall Metrics Summary',
                               data['overall'] if not data['overall'].empty else data['accuracy'],
                               window=windows[0] if windows else None)

        # Metrics by modification type (as figure instead of table)
        for window in windows:
            plot_metrics_by_modification_type(data, pdf, window)

        # Per-tool detailed pages
        add_per_tool_pages(data, pdf)

        # Optimal thresholds
        add_optimal_thresholds_page(data, pdf)

        # Score distributions
        add_score_distributions_page(data, pdf)

        # Data sources page
        add_data_sources_page(pdf, benchmark_dir)

    return output_path


def main():
    if 'snakemake' in globals():
        benchmark_dir = snakemake.params.benchmark_dir
        output_pdf = snakemake.output.pdf
    else:
        parser = argparse.ArgumentParser(description='Generate PDF benchmark report')
        parser.add_argument('--benchmark-dir', required=True,
                           help='Directory containing benchmark TSV files')
        parser.add_argument('--output', required=True,
                           help='Output PDF file path')
        args = parser.parse_args()

        benchmark_dir = args.benchmark_dir
        output_pdf = args.output

    # Load data
    print(f"Loading benchmark data from {benchmark_dir}")
    data = load_benchmark_data(benchmark_dir)

    # Generate report
    print(f"Generating PDF report: {output_pdf}")
    generate_pdf_report(data, output_pdf, benchmark_dir)
    print(f"Done! Report saved to {output_pdf}")


if __name__ == '__main__':
    main()
