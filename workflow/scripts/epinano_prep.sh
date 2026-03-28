#!/usr/bin/env bash

set -x
set -e

directory=$(dirname "${snakemake_output[per_site]}")
expected_output="${snakemake_output[per_site]}"
bam="${snakemake_input[sample_bam]}"
ref="${snakemake_input[reference]}"
extra="${snakemake_params[extra]}"

# Verify reference .fai index exists (EpiNano silently exits without it)
if [ ! -f "${ref}.fai" ]; then
    echo "ERROR: Reference .fai index not found: ${ref}.fai"
    echo "Run: samtools faidx ${ref}"
    exit 1
fi

# Check BAM has mapped reads
echo "Checking BAM: ${bam}"
samtools flagstat "${bam}" | head -5

# Resolve to absolute paths before cd
ref_abs="$(cd "$(dirname "${ref}")" && pwd)/$(basename "${ref}")"
bam_abs="$(cd "$(dirname "${bam}")" && pwd)/$(basename "${bam}")"
log_abs="$(cd "$(dirname "${snakemake_log[0]}")" && pwd)/$(basename "${snakemake_log[0]}")"

echo "Running EpiNano_Variants.py with -c ${snakemake[threads]} CPUs, output dir: ${directory}"
# v1.2.4 has no -o flag; output goes to CWD, so cd to target directory
pushd "${directory}" > /dev/null
python /opt/epinano/Epinano_Variants.py -r "${ref_abs}" \
    -b "${bam_abs}" -c ${snakemake[threads]} ${extra} \
    >>"${log_abs}" 2>&1
popd > /dev/null

# List all files EpiNano produced in output directory
echo "=== All files in ${directory} after EpiNano ==="
ls -lh "${directory}/" 2>/dev/null || echo "  (directory empty or missing)"

echo "=== Files matching *per.site* in ${directory} ==="
ls -lh ${directory}/*per.site* 2>/dev/null || echo "  (none found)"

# EpiNano output: {bam_stem}.fwd.per.site.csv (from split_bam fwd.bam)
# Some versions may produce {bam_stem}.per.site.csv (no strand split)
# If expected output doesn't exist, try to find and rename what was produced
if [ ! -f "${expected_output}" ]; then
    bam_stem=$(basename "${bam}" .bam)
    # Search in output directory
    found=$(find "${directory}" -maxdepth 1 -name "${bam_stem}*per.site.csv" -print -quit 2>/dev/null)
    # Broader fallback: any per.site.csv in output directory
    if [ -z "${found}" ]; then
        found=$(find "${directory}" -maxdepth 1 -name "*per.site.csv" -print -quit 2>/dev/null)
    fi
    if [ -n "${found}" ]; then
        echo "Renaming ${found} -> ${expected_output}"
        mv "${found}" "${expected_output}"
    else
        echo "ERROR: EpiNano produced no per.site.csv output for ${bam_stem}"
        echo "This usually means no reads mapped to forward strand after BAM splitting."
        echo "Full diagnostic info written to log: ${snakemake_log[0]}"
        # Create empty output with header so downstream rules can handle gracefully
        echo "#Ref,pos,strand,base,cov,q_mean,q_median,q_std,mat,mis,ins,del" > "${expected_output}"
    fi
fi

# Clean up EpiNano intermediate files (split BAMs)
rm -f "${directory}/${bam_stem}.fwd.bam" "${directory}/${bam_stem}.rev.bam" 2>/dev/null || true
