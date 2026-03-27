#!/usr/bin/env bash

set -x
set -e

directory=$(dirname "${snakemake_output[per_site]}")
expected_output="${snakemake_output[per_site]}"
bam="${snakemake_input[sample_bam]}"
ref="${snakemake_input[reference]}"

# Verify reference .fai index exists (EpiNano silently exits without it)
if [ ! -f "${ref}.fai" ]; then
    echo "ERROR: Reference .fai index not found: ${ref}.fai"
    echo "Run: samtools faidx ${ref}"
    exit 1
fi

# Check BAM has mapped reads
echo "Checking BAM: ${bam}"
samtools flagstat "${bam}" | head -5

python /opt/epinano/Epinano_Variants.py -r "${ref}" \
    -b "${bam}" -c ${snakemake[threads]} -o "${directory}" 2>"${snakemake_log[0]}"

# List what EpiNano actually produced
echo "Files matching *per.site* in ${directory}:"
ls -lh ${directory}/*per.site* 2>/dev/null || echo "  (none found)"

# EpiNano v1.2.5 output: {bam_stem}.fwd.per.site.csv (from split_bam fwd.bam)
# If expected output doesn't exist, try to find and rename what was produced
if [ ! -f "${expected_output}" ]; then
    # Try any per.site.csv file containing the sample name
    bam_stem=$(basename "${bam}" .bam)
    found=$(find "${directory}" -maxdepth 1 -name "${bam_stem}*per.site.csv" -print -quit 2>/dev/null)
    if [ -n "${found}" ]; then
        echo "Renaming ${found} -> ${expected_output}"
        mv "${found}" "${expected_output}"
    else
        echo "ERROR: EpiNano produced no per.site.csv output for ${bam_stem}"
        echo "This usually means no reads mapped to forward strand after BAM splitting."
        echo "Check the log: ${snakemake_log[0]}"
        # Create empty output with header so downstream rules can handle gracefully
        echo "#Ref,pos,strand,base,cov,q_mean,q_median,q_std,mat,mis,ins,del" > "${expected_output}"
    fi
fi
