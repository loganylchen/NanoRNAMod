"""
Generate visualization plots for benchmarking results.

Generates:
- Precision-Recall curves
- ROC curves
- F1 score by window bar plots
- Resource comparison plots
- Optimal threshold analysis
"""

import os
import argparse
import pandas as pd
import numpy as np
from pathlib import Path

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
except ImportError:
    print("Warning: matplotlib/seaborn not available. Install with: pip install matplotlib seaborn")
    raise


def load_accuracy_summary(path):
    """Load accuracy_summary.tsv or accuracy_summary_overall.tsv."""
    df = pd.read_csv(path, sep='\t', low_memory=False)
    if 'window' not in df.columns:
        df['window'] = 0
    return df


def load_resource_summary(path):
    """Load resource_summary.tsv."""
    return pd.read_csv(path, sep='\t', low_memory=False)


def plot_pr_curve(df, output_path, window=0):
    """Generate Precision-Recall curve comparison."""
    df_w = df[df['window'] == window].copy()

    if 'auprc' not in df.columns or df_w.empty:
        return None

    # Filter to tools with AUPRC values
    df_w = df_w.dropna(subset=['auprc'])

    if df_w.empty:
        return None

    fig, ax = plt.subplots(figsize=(10, 6))

    # Group by tool (and modification_type if present)
    if 'modification_type' in df.columns:
        # Per-modification plot
        mod_types = df_w['modification_type'].unique()
        tools = df_w['tool'].unique()

        x = np.arange(len(tools))
        width = 0.8 / len(mod_types)

        for i, mod in enumerate(mod_types):
            subset = df_w[df_w['modification_type'] == mod]
            values = [subset[subset['tool'] == t]['auprc'].values[0]
                      if t in subset['tool'].values else 0
                      for t in tools]
            ax.bar(x + i * width, values, width, label=mod)

        ax.set_xlabel('Tool')
        ax.set_ylabel('AUPRC')
        ax.set_title(f'AUPRC by Tool and Modification Type (window={window})')
        ax.set_xticks(x + width * (len(mod_types) - 1) / 2)
        ax.set_xticklabels(tools, rotation=45, ha='right')
        ax.legend()
    else:
        # Overall comparison
        tools = df_w['tool'].unique()
        auprcs = [df_w[df_w['tool'] == t]['auprc'].values[0]
                  for t in tools]

        ax.bar(tools, auprcs)
        ax.set_xlabel('Tool')
        ax.set_ylabel('AUPRC')
        ax.set_title(f'AUPRC by Tool (window={window})')
        ax.tick_params(axis='x', rotation=45)

    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def plot_roc_curve(df, output_path, window=0):
    """Generate ROC curve comparison (AUROC bar plot)."""
    df_w = df[df['window'] == window].copy()

    if 'auroc' not in df.columns or df_w.empty:
        return None

    df_w = df_w.dropna(subset=['auroc'])

    if df_w.empty:
        return None

    fig, ax = plt.subplots(figsize=(10, 6))

    if 'modification_type' in df.columns:
        mod_types = df_w['modification_type'].unique()
        tools = df_w['tool'].unique()

        x = np.arange(len(tools))
        width = 0.8 / len(mod_types)

        for i, mod in enumerate(mod_types):
            subset = df_w[df_w['modification_type'] == mod]
            values = [subset[subset['tool'] == t]['auroc'].values[0]
                      if t in subset['tool'].values else 0
                      for t in tools]
            ax.bar(x + i * width, values, width, label=mod)

        ax.set_xlabel('Tool')
        ax.set_ylabel('AUROC')
        ax.set_title(f'AUROC by Tool and Modification Type (window={window})')
        ax.set_xticks(x + width * (len(mod_types) - 1) / 2)
        ax.set_xticklabels(tools, rotation=45, ha='right')
        ax.legend()
    else:
        tools = df_w['tool'].unique()
        aurocs = [df_w[df_w['tool'] == t]['auroc'].values[0]
                  for t in tools]

        ax.bar(tools, aurocs)
        ax.set_xlabel('Tool')
        ax.set_ylabel('AUROC')
        ax.set_title(f'AUROC by Tool (window={window})')
        ax.tick_params(axis='x', rotation=45)

    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(0, 1.0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def plot_f1_by_window(df, output_path):
    """Generate F1 score comparison across windows."""
    if 'f1' not in df.columns:
        return None

    fig, ax = plt.subplots(figsize=(10, 6))

    tools = df['tool'].unique()

    if 'modification_type' in df.columns:
        # Aggregate by modification type
        mod_types = df['modification_type'].unique()
        windows = sorted(df['window'].unique())

        x = np.arange(len(windows))
        width = 0.8 / len(mod_types)

        for i, mod in enumerate(mod_types):
            values = [df[(df['modification_type'] == mod) & (df['window'] == w)]['f1'].mean()
                      for w in windows]
            ax.bar(x + i * width, values, width, label=mod)

        ax.set_xlabel('Window (nt)')
        ax.set_ylabel('F1 Score')
        ax.set_title('F1 Score by Window and Modification Type')
        ax.set_xticks(x + width * (len(mod_types) - 1) / 2)
        ax.set_xticklabels([f'{w}' for w in windows])
        ax.legend()
    else:
        windows = sorted(df['window'].unique())
        tools = df['tool'].unique()

        x = np.arange(len(windows))
        width = 0.8 / len(tools)

        for i, tool in enumerate(tools):
            values = [df[(df['tool'] == tool) & (df['window'] == w)]['f1'].values[0]
                      if len(df[(df['tool'] == tool) & (df['window'] == w)]) > 0 else 0
                      for w in windows]
            ax.bar(x + i * width, values, width, label=tool)

        ax.set_xlabel('Window (nt)')
        ax.set_ylabel('F1 Score')
        ax.set_title('F1 Score by Window')
        ax.set_xticks(x + width * (len(tools) - 1) / 2)
        ax.set_xticklabels([f'{w}' for w in windows])
        ax.legend()

    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(0, 1.0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def plot_resource_comparison(df, output_path):
    """Generate resource usage comparison plots."""
    if df.empty or 'tool' not in df.columns:
        return None

    # Build agg spec dynamically — only include columns that exist
    agg_spec = {}
    for col in ['s', 'cpu_time', 'max_rss']:
        if col in df.columns:
            agg_spec[col] = 'mean'
    for col in ['io_in', 'io_out']:
        if col in df.columns:
            agg_spec[col] = 'sum'

    if not agg_spec:
        return None

    agg = df.groupby('tool').agg(agg_spec).reset_index()

    # Determine which subplots to draw
    panels = []
    if 's' in agg.columns:
        panels.append(('s', 'Wall-Clock Time by Tool', 'Time (seconds)'))
    if 'max_rss' in agg.columns:
        panels.append(('max_rss', 'Peak Memory Usage by Tool', 'Memory (MB)'))
    if 'cpu_time' in agg.columns:
        panels.append(('cpu_time', 'Total CPU Time by Tool', 'CPU Time (seconds)'))
    if 'io_in' in agg.columns and 'io_out' in agg.columns:
        agg['total_io'] = agg['io_in'] + agg['io_out']
        panels.append(('total_io', 'Total I/O by Tool', 'I/O (MB)'))

    if not panels:
        return None

    ncols = 2
    nrows = (len(panels) + 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 5 * nrows))
    axes_flat = np.array(axes).flatten()

    for i, (col, title, xlabel) in enumerate(panels):
        ax = axes_flat[i]
        sorted_agg = agg.sort_values(col, ascending=False)
        ax.barh(sorted_agg['tool'], sorted_agg[col])
        ax.set_xlabel(xlabel)
        ax.set_title(title)
        ax.grid(axis='x', alpha=0.3)

    # Hide unused subplots
    for j in range(len(panels), len(axes_flat)):
        axes_flat[j].axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def plot_metrics_heatmap(df, output_path, window=0):
    """Generate heatmap of key metrics by tool."""
    df_w = df[df['window'] == window].copy()

    if df_w.empty:
        return None

    # Select key metrics
    metrics = ['precision', 'recall', 'f1']
    if 'auprc' in df_w.columns:
        metrics.append('auprc')
    if 'auroc' in df_w.columns:
        metrics.append('auroc')

    # Filter to columns that exist
    metrics = [m for m in metrics if m in df_w.columns]

    if 'modification_type' in df.columns:
        # Pivot table for heatmap
        pivot = df_w.pivot_table(
            index='tool',
            columns='modification_type',
            values=metrics[0]  # Use first metric for simplicity
        )
    else:
        pivot = df_w.set_index('tool')[metrics]

    if pivot.empty:
        return None

    fig, ax = plt.subplots(figsize=(10, 6))

    if 'modification_type' in df.columns:
        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn',
                   vmin=0, vmax=1, ax=ax, cbar_kws={'label': metrics[0].title()})
        ax.set_title(f'{metrics[0].title()} by Tool and Modification Type (window={window})')
        ax.set_ylabel('Tool')
        ax.set_xlabel('Modification Type')
    else:
        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn',
                   vmin=0, vmax=1, ax=ax, cbar_kws={'label': 'Score'})
        ax.set_title(f'Metrics by Tool (window={window})')
        ax.set_ylabel('Tool')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def plot_optimal_score_columns(optimal_df, output_path):
    """Visualize optimal score columns selected for each tool.

    Shows which score column is best for each tool and the F1 achieved.
    """
    if optimal_df.empty:
        return None

    required_cols = ['tool', 'score_column', 'f1']
    if not all(c in optimal_df.columns for c in required_cols):
        return None

    fig, ax = plt.subplots(figsize=(12, max(6, len(optimal_df) * 0.5)))

    # Sort by F1 score
    df = optimal_df.sort_values('f1', ascending=True)

    # Create horizontal bar chart
    y_pos = np.arange(len(df))
    bars = ax.barh(y_pos, df['f1'], color='steelblue', edgecolor='navy')

    # Add score column names as annotations
    for i, (idx, row) in enumerate(df.iterrows()):
        ax.text(row['f1'] + 0.01, i,
                f"{row['score_column']} (F1={row['f1']:.3f})",
                va='center', fontsize=9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df['tool'])
    ax.set_xlabel('F1 Score')
    ax.set_title('Optimal Score Column Selection by Tool')
    ax.set_xlim(0, 1.1)
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def plot_optimal_thresholds(threshold_df, output_path):
    """Visualize optimal thresholds for each tool.

    Shows the optimal threshold (in -log10 scale) and the original p-value.
    """
    if threshold_df.empty:
        return None

    required_cols = ['tool', 'optimal_threshold']
    if not all(c in threshold_df.columns for c in required_cols):
        return None

    # Filter to valid thresholds
    df = threshold_df.dropna(subset=['optimal_threshold']).copy()
    if df.empty:
        return None

    # Sort by threshold
    df = df.sort_values('optimal_threshold', ascending=True)

    fig, ax = plt.subplots(figsize=(12, max(6, len(df) * 0.5)))

    y_pos = np.arange(len(df))

    # Plot thresholds
    bars = ax.barh(y_pos, df['optimal_threshold'],
                   color=['green' if t > 1.3 else 'orange' if t > 0.5 else 'red'
                          for t in df['optimal_threshold']],
                   edgecolor='navy')

    # Add annotations with original threshold if available
    for i, (idx, row) in enumerate(df.iterrows()):
        threshold_text = f"{row['optimal_threshold']:.2f}"
        if 'original_threshold' in df.columns and pd.notna(row.get('original_threshold')):
            threshold_text += f" (p={row['original_threshold']:.2e})"
        ax.text(row['optimal_threshold'] + 0.05, i, threshold_text,
                va='center', fontsize=9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df['tool'])
    ax.set_xlabel('Optimal Threshold (-log10 p-value)')
    ax.set_title('Optimal Score Thresholds by Tool')
    ax.axvline(x=1.3, color='gray', linestyle='--', alpha=0.5, label='p=0.05')
    ax.axvline(x=2.0, color='gray', linestyle=':', alpha=0.5, label='p=0.01')
    ax.legend(loc='lower right', fontsize=8)
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def plot_threshold_performance_curves(threshold_eval_df, output_path):
    """Plot precision, recall, F1 as function of threshold for each tool."""
    if threshold_eval_df.empty:
        return None

    required_cols = ['tool', 'threshold', 'precision', 'recall', 'f1']
    if not all(c in threshold_eval_df.columns for c in required_cols):
        return None

    tools = threshold_eval_df['tool'].unique()
    if len(tools) == 0:
        return None

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    metrics = ['precision', 'recall', 'f1']
    colors = plt.cm.tab10(np.linspace(0, 1, len(tools)))
    color_map = dict(zip(tools, colors))

    for ax, metric in zip(axes, metrics):
        for tool in tools:
            tool_df = threshold_eval_df[threshold_eval_df['tool'] == tool].sort_values('threshold')
            if not tool_df.empty:
                ax.plot(tool_df['threshold'], tool_df[metric],
                        label=tool, color=color_map[tool], alpha=0.7, linewidth=2)

        ax.set_xlabel('Threshold (-log10 p-value)')
        ax.set_ylabel(metric.capitalize())
        ax.set_title(f'{metric.capitalize()} vs Threshold')
        ax.grid(alpha=0.3)
        ax.set_ylim(0, 1)

    # Add legend to last plot
    axes[-1].legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def plot_score_distributions(score_dist_df, output_path):
    """Visualize score distribution statistics for each tool."""
    if score_dist_df.empty:
        return None

    required_cols = ['tool', 'mean', 'std', 'min', 'max']
    if not all(c in score_dist_df.columns for c in required_cols):
        return None

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    df = score_dist_df.sort_values('mean', ascending=False)

    # Left plot: Mean ± std
    ax = axes[0]
    y_pos = np.arange(len(df))
    ax.barh(y_pos, df['mean'], xerr=df['std'], color='steelblue',
            edgecolor='navy', alpha=0.7, capsize=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df['tool'])
    ax.set_xlabel('Score Value')
    ax.set_title('Score Distribution: Mean ± Std')
    ax.grid(axis='x', alpha=0.3)

    # Right plot: Range (min to max)
    ax = axes[1]
    for i, (idx, row) in enumerate(df.iterrows()):
        ax.plot([row['min'], row['max']], [i, i], 'o-', color='steelblue', alpha=0.7)
        ax.plot(row['mean'], i, 'ro', markersize=8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df['tool'])
    ax.set_xlabel('Score Value')
    ax.set_title('Score Range (min-mean-max)')
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def plot_score_column_comparison(all_scores_df, output_path):
    """Compare performance across different score columns for each tool."""
    if all_scores_df.empty:
        return None

    required_cols = ['tool', 'score_column', 'f1']
    if not all(c in all_scores_df.columns for c in required_cols):
        return None

    tools = all_scores_df['tool'].unique()
    n_tools = len(tools)

    if n_tools == 0:
        return None

    fig, axes = plt.subplots(n_tools, 1, figsize=(10, 4 * n_tools), squeeze=False)
    axes = axes.flatten()

    for i, tool in enumerate(tools):
        tool_df = all_scores_df[all_scores_df['tool'] == tool].sort_values('f1', ascending=True)

        if tool_df.empty:
            continue

        ax = axes[i]
        y_pos = np.arange(len(tool_df))

        # Color by whether it's the best
        best_f1 = tool_df['f1'].max()
        colors = ['green' if f == best_f1 else 'steelblue' for f in tool_df['f1']]

        ax.barh(y_pos, tool_df['f1'], color=colors, edgecolor='navy')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(tool_df['score_column'])
        ax.set_xlabel('F1 Score')
        ax.set_title(f'{tool}: Score Column Comparison')
        ax.set_xlim(0, 1)
        ax.grid(axis='x', alpha=0.3)

        # Mark best
        best_idx = tool_df['f1'].idxmax()
        best_pos = y_pos[tool_df.index.get_loc(best_idx)]
        ax.text(1.02, best_pos, '★ Best', va='center', fontsize=9, color='green')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def generate_html_report(accuracy_df, resource_df, figures_dir, output_path):
    """Generate self-contained HTML benchmark report."""

    def safe_float(val):
        """Convert value to float, return NaN if invalid."""
        try:
            return float(val) if not pd.isna(val) else None
        except (TypeError, ValueError):
            return None

    # Prepare accuracy data
    accuracy_rows = []
    for _, row in accuracy_df.iterrows():
        acc_row = {
            'tool': row.get('tool', ''),
            'window': int(row.get('window', 0)),
            'precision': safe_float(row.get('precision')),
            'recall': safe_float(row.get('recall')),
            'f1': safe_float(row.get('f1')),
            'tp': int(row.get('tp', 0)) if pd.notna(row.get('tp')) else 0,
            'fp': int(row.get('fp', 0)) if pd.notna(row.get('fp')) else 0,
            'fn': int(row.get('fn', 0)) if pd.notna(row.get('fn')) else 0,
            'tn': int(row.get('tn', 0)) if pd.notna(row.get('tn')) else 0,
            'auprc': safe_float(row.get('auprc')),
            'auroc': safe_float(row.get('auroc')),
            'mcc': safe_float(row.get('mcc')),
        }
        if 'modification_type' in accuracy_df.columns:
            acc_row['modification_type'] = row.get('modification_type', '')
        accuracy_rows.append(acc_row)

    # Prepare resource data
    resource_rows = []
    for _, row in resource_df.iterrows():
        resource_rows.append({
            'tool': row.get('tool', ''),
            'sample': row.get('sample', ''),
            'time_s': safe_float(row.get('s')),
            'max_rss_mb': safe_float(row.get('max_rss')),
            'cpu_time_s': safe_float(row.get('cpu_time')),
            'io_in_mb': safe_float(row.get('io_in')),
            'io_out_mb': safe_float(row.get('io_out')),
        })

    # Get list of generated figures
    figures = []
    if os.path.exists(figures_dir):
        for f in os.listdir(figures_dir):
            if f.endswith(('.svg', '.png')):
                figures.append(f)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NanoRNAMod Benchmark Report</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
               margin: 0; padding: 20px; background: #f5f5f5; }}
        .container {{ max-width: 1400px; margin: 0 auto; background: white; padding: 30px;
                     border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background: #3498db; color: white; font-weight: 600; }}
        tr:hover {{ background: #f8f9fa; }}
        .metric {{ display: inline-block; padding: 6px 12px; margin: 4px;
                  border-radius: 4px; font-weight: 600; }}
        .metric-high {{ background: #27ae60; color: white; }}
        .metric-med {{ background: #f39c12; color: white; }}
        .metric-low {{ background: #e74c3c; color: white; }}
        .figure {{ text-align: center; margin: 30px 0; }}
        .figure img {{ max-width: 100%; border: 1px solid #ddd; border-radius: 4px; }}
        .tab-buttons {{ margin: 20px 0; }}
        .tab-btn {{ padding: 10px 20px; margin-right: 5px; border: none;
                   background: #ecf0f1; cursor: pointer; border-radius: 4px; }}
        .tab-btn.active {{ background: #3498db; color: white; }}
        .tab-content {{ display: none; }}
        .tab-content.active {{ display: block; }}
        .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                 gap: 15px; margin: 20px 0; }}
        .card {{ background: #ecf0f1; padding: 15px; border-radius: 8px; }}
        .card-label {{ font-size: 12px; color: #7f8c8d; }}
        .card-value {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 NanoRNAMod Benchmark Report</h1>

        <div class="tab-buttons">
            <button class="tab-btn active" onclick="showTab('overview')">Overview</button>
            <button class="tab-btn" onclick="showTab('accuracy')">Accuracy</button>
            <button class="tab-btn" onclick="showTab('resources')">Resources</button>
            <button class="tab-btn" onclick="showTab('figures')">Figures</button>
        </div>

        <div id="overview" class="tab-content active">
            <h2>Summary</h2>
            <div class="grid">
                <div class="card">
                    <div class="card-label">Tools Evaluated</div>
                    <div class="card-value">{len(set(r['tool'] for r in accuracy_rows))}</div>
                </div>
                <div class="card">
                    <div class="card-label">Total Predictions</div>
                    <div class="card-value">{sum(r.get('called_sites', 0) for r in accuracy_rows)}</div>
                </div>
                <div class="card">
                    <div class="card-label">Windows Tested</div>
                    <div class="card-value">{len(set(r['window'] for r in accuracy_rows))}</div>
                </div>
            </div>
        </div>

        <div id="accuracy" class="tab-content">
            <h2>Accuracy Metrics</h2>
            <table id="accuracyTable">
                <thead>
                    <tr>
                        <th>Tool</th>
"""

    if 'modification_type' in accuracy_df.columns:
        html += "                        <th>Modification</th>"

    html += """                        <th>Window</th>
                        <th>Precision</th>
                        <th>Recall</th>
                        <th>F1</th>
                        <th>TP</th>
                        <th>FP</th>
                        <th>FN</th>
                        <th>AUROC</th>
                        <th>AUPRC</th>
                    </tr>
                </thead>
                <tbody>
"""

    for row in accuracy_rows:
        f1 = row['f1'] or 0
        f1_class = 'metric-high' if f1 >= 0.8 else 'metric-med' if f1 >= 0.5 else 'metric-low'

        html += f"""                    <tr>
                        <td>{row['tool']}</td>
"""
        if 'modification_type' in row:
            html += f"                        <td>{row.get('modification_type', '')}</td>\n"

        # Format metrics safely (handle None/NaN values)
        prec_str = f"{row['precision']:.3f}" if pd.notna(row.get('precision')) else 'N/A'
        recall_str = f"{row['recall']:.3f}" if pd.notna(row.get('recall')) else 'N/A'
        f1_str = f"{row['f1']:.3f}" if pd.notna(row.get('f1')) else 'N/A'
        auroc_str = f"{row['auroc']:.3f}" if pd.notna(row.get('auroc')) else 'N/A'
        auprc_str = f"{row['auprc']:.3f}" if pd.notna(row.get('auprc')) else 'N/A'

        html += f"""                        <td>{row['window']}</td>
                        <td><span class="metric {f1_class}">{prec_str}</span></td>
                        <td><span class="metric {f1_class}">{recall_str}</span></td>
                        <td><span class="metric {f1_class}">{f1_str}</span></td>
                        <td>{row['tp']}</td>
                        <td>{row['fp']}</td>
                        <td>{row['fn']}</td>
                        <td>{auroc_str}</td>
                        <td>{auprc_str}</td>
                    </tr>
"""

    html += """                </tbody>
            </table>
        </div>

        <div id="resources" class="tab-content">
            <h2>Resource Usage</h2>
            <table>
                <thead>
                    <tr>
                        <th>Tool</th>
                        <th>Sample</th>
                        <th>Time (s)</th>
                        <th>Memory (MB)</th>
                        <th>CPU Time (s)</th>
                        <th>I/O In (MB)</th>
                        <th>I/O Out (MB)</th>
                    </tr>
                </thead>
                <tbody>
"""

    for row in resource_rows:
        # Format metrics safely
        time_str = f"{row['time_s']:.2f}" if pd.notna(row.get('time_s')) else 'N/A'
        mem_str = f"{row['max_rss_mb']:.0f}" if pd.notna(row.get('max_rss_mb')) else 'N/A'
        cpu_str = f"{row['cpu_time_s']:.2f}" if pd.notna(row.get('cpu_time_s')) else 'N/A'
        io_in_str = f"{row['io_in_mb']:.2f}" if pd.notna(row.get('io_in_mb')) else 'N/A'
        io_out_str = f"{row['io_out_mb']:.2f}" if pd.notna(row.get('io_out_mb')) else 'N/A'

        html += f"""                    <tr>
                        <td>{row['tool']}</td>
                        <td>{row['sample']}</td>
                        <td>{time_str}</td>
                        <td>{mem_str}</td>
                        <td>{cpu_str}</td>
                        <td>{io_in_str}</td>
                        <td>{io_out_str}</td>
                    </tr>
"""

    html += """                </tbody>
            </table>
        </div>

        <div id="figures" class="tab-content">
            <h2>Visualization Figures</h2>
"""

    for fig in figures:
        html += f"""            <div class="figure">
                <img src="figures/{fig}" alt="{fig}">
                <p>{fig.replace('_', ' ').replace('.svg', '').title()}</p>
            </div>
"""

    html += """        </div>
    </div>

    <script>
        function showTab(tabId) {
            document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
            document.querySelectorAll('.tab-btn').forEach(el => el.classList.remove('active'));
            document.getElementById(tabId).classList.add('active');
            event.target.classList.add('active');
        }
    </script>
</body>
</html>
"""

    with open(output_path, 'w') as f:
        f.write(html)

    return output_path


def main():
    import sys

    # When called by Snakemake
    if 'snakemake' in globals():
        accuracy_file = snakemake.input.accuracy
        resource_file = snakemake.input.resource
        # Get output directory from params (parent of html output)
        output_dir = snakemake.params.output_dir
        html_output = snakemake.output.html
        figures_dir = os.path.join(output_dir, 'figures')
        os.makedirs(figures_dir, exist_ok=True)

        accuracy_df = load_accuracy_summary(accuracy_file)
        resource_df = load_resource_summary(resource_file)

        # Generate plots
        plots = []
        windows = sorted(accuracy_df['window'].unique()) if 'window' in accuracy_df.columns else [0]

        for window in windows:
            pr_path = os.path.join(figures_dir, f'pr_curve_window{window}.svg')
            if plot_pr_curve(accuracy_df, pr_path, window):
                plots.append(pr_path)

            roc_path = os.path.join(figures_dir, f'roc_curve_window{window}.svg')
            if plot_roc_curve(accuracy_df, roc_path, window):
                plots.append(roc_path)

            heatmap_path = os.path.join(figures_dir, f'heatmap_window{window}.svg')
            if plot_metrics_heatmap(accuracy_df, heatmap_path, window):
                plots.append(heatmap_path)

        f1_path = os.path.join(figures_dir, 'f1_by_window.svg')
        if plot_f1_by_window(accuracy_df, f1_path):
            plots.append(f1_path)

        resource_path = os.path.join(figures_dir, 'resource_comparison.svg')
        if plot_resource_comparison(resource_df, resource_path):
            plots.append(resource_path)

        # Generate optimal parameter visualizations if data exists
        benchmark_dir = os.path.dirname(output_dir) if os.path.basename(output_dir) == 'viz' else output_dir

        # Optimal score columns
        optimal_scores_path = os.path.join(benchmark_dir, 'optimal_score_per_tool.tsv')
        if os.path.exists(optimal_scores_path):
            optimal_df = pd.read_csv(optimal_scores_path, sep='\t')
            opt_scores_plot = os.path.join(figures_dir, 'optimal_score_columns.svg')
            if plot_optimal_score_columns(optimal_df, opt_scores_plot):
                plots.append(opt_scores_plot)

        # All scores evaluation (for comparison plots)
        all_scores_path = os.path.join(benchmark_dir, 'all_scores_evaluation.tsv')
        if os.path.exists(all_scores_path):
            all_scores_df = pd.read_csv(all_scores_path, sep='\t')
            scores_comp_plot = os.path.join(figures_dir, 'score_column_comparison.svg')
            if plot_score_column_comparison(all_scores_df, scores_comp_plot):
                plots.append(scores_comp_plot)

        # Optimal thresholds
        optimal_thresh_path = os.path.join(benchmark_dir, 'optimal_thresholds_detailed.tsv')
        if os.path.exists(optimal_thresh_path):
            thresh_df = pd.read_csv(optimal_thresh_path, sep='\t')
            thresh_plot = os.path.join(figures_dir, 'optimal_thresholds.svg')
            if plot_optimal_thresholds(thresh_df, thresh_plot):
                plots.append(thresh_plot)

        # Generate HTML report
        generate_html_report(accuracy_df, resource_df, figures_dir, html_output)

        print(f"Generated {len(plots)} plots and HTML report")

    # When called from command line
    else:
        parser = argparse.ArgumentParser(description='Generate benchmark visualizations')
        parser.add_argument('--accuracy', required=True, help='accuracy_summary.tsv path')
        parser.add_argument('--resource', required=True, help='resource_summary.tsv path')
        parser.add_argument('--output', required=True, help='Output directory')
        parser.add_argument('--benchmark-dir', default=None, help='Benchmark directory for optimal params')
        args = parser.parse_args()

        os.makedirs(os.path.join(args.output, 'figures'), exist_ok=True)

        accuracy_df = load_accuracy_summary(args.accuracy)
        resource_df = load_resource_summary(args.resource)

        figures_dir = os.path.join(args.output, 'figures')

        for window in sorted(accuracy_df['window'].unique()) if 'window' in accuracy_df.columns else [0]:
            plot_pr_curve(accuracy_df, os.path.join(figures_dir, f'pr_curve_window{window}.svg'), window)
            plot_roc_curve(accuracy_df, os.path.join(figures_dir, f'roc_curve_window{window}.svg'), window)
            plot_metrics_heatmap(accuracy_df, os.path.join(figures_dir, f'heatmap_window{window}.svg'), window)

        plot_f1_by_window(accuracy_df, os.path.join(figures_dir, 'f1_by_window.svg'))
        plot_resource_comparison(resource_df, os.path.join(figures_dir, 'resource_comparison.svg'))

        # Generate optimal parameter visualizations
        benchmark_dir = args.benchmark_dir or args.output
        optimal_scores_path = os.path.join(benchmark_dir, 'optimal_score_per_tool.tsv')
        if os.path.exists(optimal_scores_path):
            optimal_df = pd.read_csv(optimal_scores_path, sep='\t')
            plot_optimal_score_columns(optimal_df, os.path.join(figures_dir, 'optimal_score_columns.svg'))

        all_scores_path = os.path.join(benchmark_dir, 'all_scores_evaluation.tsv')
        if os.path.exists(all_scores_path):
            all_scores_df = pd.read_csv(all_scores_path, sep='\t')
            plot_score_column_comparison(all_scores_df, os.path.join(figures_dir, 'score_column_comparison.svg'))

        optimal_thresh_path = os.path.join(benchmark_dir, 'optimal_thresholds_detailed.tsv')
        if os.path.exists(optimal_thresh_path):
            thresh_df = pd.read_csv(optimal_thresh_path, sep='\t')
            plot_optimal_thresholds(thresh_df, os.path.join(figures_dir, 'optimal_thresholds.svg'))

        generate_html_report(accuracy_df, resource_df, figures_dir,
                            os.path.join(args.output, 'benchmark_report.html'))

        print(f"Generated plots in {figures_dir}")
        print(f"Generated HTML report at {os.path.join(args.output, 'benchmark_report.html')}")


if __name__ == '__main__':
    main()
