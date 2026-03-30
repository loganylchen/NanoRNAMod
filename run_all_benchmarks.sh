#!/bin/bash
# Run NanoRNAMod benchmarking across all ecoli rRNA datasets.
#
# Each dataset at ecoli_rRNAs/{depth}/{stoichiometry}/data is symlinked as ./data,
# snakemake is run, then results are renamed to ecoliRNA{depth}_{stoichiometry}.
#
# Usage: bash run_all_benchmarks.sh

set -euo pipefail

BASE_DATA="/autofs/NAS25_Shared/logan/Baleen_Stoichiometry_Testdata/ecoli_rRNAs"
PROJECT="RNAModProject"  # must match config/config.yaml project name
LOGDIR="run_logs"

mkdir -p "$LOGDIR"

# Collect all datasets (data may be a real dir or a symlink)
datasets=( $(find "$BASE_DATA" -mindepth 3 -maxdepth 3 \( -type d -o -type l \) -name data | sort) )

echo "Found ${#datasets[@]} datasets to process"
echo "================================================"

for dataset_path in "${datasets[@]}"; do
    # Extract depth and stoichiometry from path
    # e.g. /...ecoli_rRNAs/500/0.1/data -> depth=500, stoich=0.1
    parent="$(dirname "$dataset_path")"
    stoich="$(basename "$parent")"
    depth="$(basename "$(dirname "$parent")")"

    # Build output name: ecoliRNA500_01 (replace . with nothing)
    stoich_clean="${stoich//.}"
    run_name="ecoliRNA${depth}_${stoich_clean}"

    echo ""
    echo "================================================"
    echo "Processing: depth=${depth}, stoichiometry=${stoich}"
    echo "  Source: ${dataset_path}"
    echo "  Output: ${run_name}"
    echo "================================================"

    # Skip if results already exist
    if [ -d "${run_name}" ]; then
        echo "  SKIP: ${run_name} already exists"
        continue
    fi

    # Remove old data symlink if present
    if [ -L "data" ]; then
        rm data
    elif [ -e "data" ]; then
        echo "  ERROR: ./data exists but is not a symlink. Skipping."
        continue
    fi

    # Create symlink
    ln -s "$dataset_path" data
    echo "  Symlinked: data -> ${dataset_path}"

    # Clean previous snakemake state for a fresh run
    rm -rf .snakemake/incomplete .snakemake/locks

    # Run snakemake
    echo "  Running snakemake..."
    snakemake \
        --jobs 15 -k \
        --slurm \
        --jobname "{rule}" \
        --default-resources slurm_account=users mem_mb=180000 \
        --resources gpu=1 \
        --cores 10 \
        --rerun-incomplete \
        --show-failed-logs \
        -p \
        --use-singularity \
        --singularity-args '-B /autofs/NAS25_Shared --nv' \
        2>&1 | tee "${LOGDIR}/${run_name}.log"

    exit_code=${PIPESTATUS[0]}

    if [ $exit_code -ne 0 ]; then
        echo "  WARNING: snakemake exited with code ${exit_code} for ${run_name}"
    fi

    # Remove data symlink
    rm data

    # Rename project results directory
    if [ -d "$PROJECT" ]; then
        mv "$PROJECT" "$run_name"
        echo "  Renamed: ${PROJECT} -> ${run_name}"
    else
        echo "  WARNING: ${PROJECT} directory not found after run"
    fi

    # Also rename logs and benchmarks if they reference the project
    if [ -d "logs/$PROJECT" ]; then
        mkdir -p "logs/${run_name}"
        mv "logs/$PROJECT"/* "logs/${run_name}/" 2>/dev/null || true
        rmdir "logs/$PROJECT" 2>/dev/null || true
    fi
    if [ -d "benchmarks/$PROJECT" ]; then
        mkdir -p "benchmarks/${run_name}"
        mv "benchmarks/$PROJECT"/* "benchmarks/${run_name}/" 2>/dev/null || true
        rmdir "benchmarks/$PROJECT" 2>/dev/null || true
    fi

    echo "  Done: ${run_name}"
done

echo ""
echo "================================================"
echo "All datasets processed."
echo "================================================"
