#!/usr/bin/env python3
"""
Reorganize benchmark outputs by negative control strategy.

This utility script reorganizes the flat benchmark output structure into
strategy-specific subdirectories. It can be run after the workflow completes
to create a more organized view of the results.

Current flat structure:
    benchmarks/
    ├── accuracy_summary.tsv
    ├── accuracy_summary_kmer_negatives.tsv
    ├── accuracy_summary_same_base_negatives.tsv
    ├── kmer_negative_sites.tsv
    ├── same_base_negative_sites.tsv
    └── ...

Reorganized structure:
    benchmarks/
    ├── shared/
    │   ├── resource_summary.tsv
    │   ├── threshold_evaluation.tsv
    │   └── viz/
    ├── default_negatives/
    │   ├── accuracy_summary.tsv
    │   └── called_sites_summary.tsv
    ├── kmer_negatives/
    │   ├── accuracy_summary.tsv
    │   └── negative_sites.tsv
    └── same_base_negatives/
        ├── accuracy_summary.tsv
        └── negative_sites.tsv

Usage:
    python organize_benchmark_outputs.py --benchmark-dir results/benchmarks [--reorganize] [--symlinks]
"""

import os
import sys
import argparse
import shutil
from pathlib import Path


# File categorization by negative control strategy
STRATEGY_FILES = {
    'shared': [
        # Resource metrics (shared across all tools)
        'resource_summary.tsv',
        'resource_by_tool.tsv',
        'resource_by_tool.pdf',
        # Multi-threshold evaluation (shared)
        'threshold_evaluation.tsv',
        'optimal_thresholds_detailed.tsv',
        'score_distributions.tsv',
        'optimal_score_per_tool.tsv',
        'all_scores_evaluation.tsv',
        # Optimal thresholds (shared)
        'optimal_thresholds.tsv',
        'detailed_predictions.tsv',
        'detailed_truth.tsv',
        # Count tables (shared across strategies)
        'called_sites_by_comparison.tsv',
        'called_sites_summary.tsv',
        # Overall summaries (aggregated)
        'accuracy_summary_overall.tsv',
        'accuracy_summary_by_comparison.tsv',
        # Visualization outputs
        'viz/',
        'figures/',
        'benchmark_report.pdf',
        'benchmark_report.html',
        # Touch files
        '.benchmark_complete',
    ],
    'default_negatives': [
        'accuracy_summary.tsv',
    ],
    'kmer_negatives': [
        'accuracy_summary_kmer_negatives.tsv',
        'kmer_negative_sites.tsv',
    ],
    'same_base_negatives': [
        'accuracy_summary_same_base_negatives.tsv',
        'same_base_negative_sites.tsv',
    ],
}


def get_strategy_for_file(filename):
    """Determine which strategy a file belongs to."""
    for strategy, files in STRATEGY_FILES.items():
        for pattern in files:
            if pattern.endswith('/') and filename.startswith(pattern.rstrip('/')):
                return strategy
            if filename == pattern:
                return strategy
    return None


def organize_benchmark_outputs(benchmark_dir, use_symlinks=False, dry_run=False):
    """
    Reorganize benchmark outputs into strategy-specific subdirectories.

    Args:
        benchmark_dir: Path to benchmarks directory
        use_symlinks: If True, create symlinks instead of moving files
        dry_run: If True, only print what would be done
    """
    benchmark_path = Path(benchmark_dir)

    if not benchmark_path.exists():
        print(f"Error: Benchmark directory does not exist: {benchmark_dir}")
        return False

    # Track what was organized
    organized = {strategy: [] for strategy in STRATEGY_FILES}
    unorganized = []

    # Get all files in benchmark directory
    all_items = list(benchmark_path.iterdir())

    for item in all_items:
        name = item.name
        strategy = get_strategy_for_file(name)

        if strategy:
            organized[strategy].append(name)
        else:
            # Check if it's already in a subdirectory
            if item.is_dir() and name in STRATEGY_FILES:
                continue
            unorganized.append(name)

    if dry_run:
        print("\n=== Dry Run - Would organize as follows ===\n")
        for strategy, files in organized.items():
            if files:
                print(f"{strategy}/")
                for f in files:
                    print(f"  {f}")
        if unorganized:
            print(f"\nUnorganized (would stay in root):")
            for f in unorganized:
                print(f"  {f}")
        return True

    # Create strategy directories
    for strategy in STRATEGY_FILES:
        strategy_dir = benchmark_path / strategy
        strategy_dir.mkdir(exist_ok=True)

    # Move/copy files to strategy directories
    for strategy, files in organized.items():
        strategy_dir = benchmark_path / strategy

        for filename in files:
            src = benchmark_path / filename
            if not src.exists():
                continue

            dst = strategy_dir / filename

            # Handle renamed files (remove strategy suffix)
            if strategy != 'shared' and strategy != 'default_negatives':
                # Rename accuracy_summary_*_negatives.tsv -> accuracy_summary.tsv
                if '_negatives' in filename:
                    new_name = filename.replace('_kmer_negatives', '').replace('_same_base_negatives', '')
                    dst = strategy_dir / new_name

            if use_symlinks:
                if dst.exists():
                    dst.unlink()
                os.symlink(os.path.relpath(src, strategy_dir), dst)
                print(f"Linked: {filename} -> {strategy}/{dst.name}")
            else:
                if dst.exists():
                    dst.unlink()
                shutil.move(str(src), str(dst))
                print(f"Moved: {filename} -> {strategy}/{dst.name}")

    print(f"\nOrganization complete!")
    print(f"  Shared files: {len(organized['shared'])}")
    print(f"  Default negatives: {len(organized['default_negatives'])}")
    print(f"  K-mer negatives: {len(organized['kmer_negatives'])}")
    print(f"  Same-base negatives: {len(organized['same_base_negatives'])}")

    if unorganized:
        print(f"\nFiles not organized (stayed in root): {len(unorganized)}")
        for f in unorganized[:10]:
            print(f"  {f}")
        if len(unorganized) > 10:
            print(f"  ... and {len(unorganized) - 10} more")

    return True


def restore_flat_structure(benchmark_dir, dry_run=False):
    """
    Restore the flat structure by moving files back to root.

    Args:
        benchmark_dir: Path to benchmarks directory
        dry_run: If True, only print what would be done
    """
    benchmark_path = Path(benchmark_dir)

    if dry_run:
        print("\n=== Dry Run - Would restore flat structure ===\n")
        for strategy in STRATEGY_FILES:
            strategy_dir = benchmark_path / strategy
            if strategy_dir.exists():
                for item in strategy_dir.iterdir():
                    print(f"Would move: {strategy}/{item.name} -> {item.name}")
        return True

    for strategy in STRATEGY_FILES:
        strategy_dir = benchmark_path / strategy

        if not strategy_dir.exists():
            continue

        for item in list(strategy_dir.iterdir()):
            if item.is_file():
                # Restore original names for strategy-specific files
                name = item.name
                if strategy == 'kmer_negatives':
                    if name == 'accuracy_summary.tsv':
                        name = 'accuracy_summary_kmer_negatives.tsv'
                    elif name == 'negative_sites.tsv':
                        name = 'kmer_negative_sites.tsv'
                elif strategy == 'same_base_negatives':
                    if name == 'accuracy_summary.tsv':
                        name = 'accuracy_summary_same_base_negatives.tsv'
                    elif name == 'negative_sites.tsv':
                        name = 'same_base_negative_sites.tsv'

                dst = benchmark_path / name
                if dst.exists():
                    dst.unlink()
                shutil.move(str(item), str(dst))
                print(f"Restored: {strategy}/{item.name} -> {name}")

        # Remove empty directory
        if strategy_dir.exists() and not list(strategy_dir.iterdir()):
            strategy_dir.rmdir()

    print("\nFlat structure restored!")
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Organize benchmark outputs by negative control strategy',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Preview reorganization
  python organize_benchmark_outputs.py --benchmark-dir results/benchmarks --dry-run

  # Reorganize using symlinks (preserves original files)
  python organize_benchmark_outputs.py --benchmark-dir results/benchmarks --symlinks

  # Reorganize by moving files
  python organize_benchmark_outputs.py --benchmark-dir results/benchmarks --reorganize

  # Restore flat structure
  python organize_benchmark_outputs.py --benchmark-dir results/benchmarks --restore
"""
    )

    parser.add_argument('--benchmark-dir', required=True,
                        help='Path to benchmarks directory')
    parser.add_argument('--dry-run', action='store_true',
                        help='Preview changes without making them')
    parser.add_argument('--symlinks', action='store_true',
                        help='Create symlinks instead of moving files')
    parser.add_argument('--reorganize', action='store_true',
                        help='Actually reorganize the files (default is dry-run)')
    parser.add_argument('--restore', action='store_true',
                        help='Restore flat structure from organized')

    args = parser.parse_args()

    if args.restore:
        return restore_flat_structure(args.benchmark_dir, dry_run=not args.reorganize)
    else:
        return organize_benchmark_outputs(
            args.benchmark_dir,
            use_symlinks=args.symlinks,
            dry_run=not args.reorganize
        )


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
