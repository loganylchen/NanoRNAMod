#!/usr/bin/env bash

set -x
set -e

expected_output="${snakemake_output[per_site]}"
bam="${snakemake_input[sample_bam]}"
ref="${snakemake_input[reference]}"
extra="${snakemake_params[extra]}"
directory=$(dirname "${expected_output}")

# Resolve to absolute paths for subshell
ref_abs="$(cd "$(dirname "${ref}")" && pwd)/$(basename "${ref}")"
bam_abs="$(cd "$(dirname "${bam}")" && pwd)/$(basename "${bam}")"

# v1.2.4 has no -o flag; output goes to CWD, so cd to target directory via subshell
# Run with explicit error capture so failures are visible
( cd "${directory}" && \
  python /opt/epinano/Epinano_Variants.py -R "${ref_abs}" \
      -s /opt/epinano/misc/sam2tsv.jar \
      -b "${bam_abs}" -t ${snakemake[threads]} ${extra} ) \
    2>&1 | tee -a "${snakemake_log[0]}"

# Check pipeline exit status (not tee's exit status)
epinano_exit=${PIPESTATUS[0]}
if [ ${epinano_exit} -ne 0 ]; then
    echo "ERROR: EpiNano_Variants.py failed with exit code ${epinano_exit}" | tee -a "${snakemake_log[0]}"
    exit ${epinano_exit}
fi

# Verify output exists
if [ ! -f "${expected_output}" ]; then
    echo "ERROR: EpiNano did not produce expected output: ${expected_output}" | tee -a "${snakemake_log[0]}"
    # List what was actually produced for debugging
    echo "Files in ${directory}:" | tee -a "${snakemake_log[0]}"
    ls -la "${directory}"/*.per.site*.csv 2>/dev/null | tee -a "${snakemake_log[0]}" || echo "  (no .per.site*.csv files found)" | tee -a "${snakemake_log[0]}"
    exit 1
fi
