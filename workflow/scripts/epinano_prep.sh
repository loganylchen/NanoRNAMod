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
( cd "${directory}" && \
  python /opt/epinano/Epinano_Variants.py -R "${ref_abs}" \
      -s /opt/epinano/misc/sam2tsv.jar \
      -b "${bam_abs}" -t ${snakemake[threads]} ${extra} ) \
    >>"${snakemake_log[0]}" 2>&1

# Verify output exists
if [ ! -f "${expected_output}" ]; then
    echo "ERROR: EpiNano did not produce expected output: ${expected_output}"
    exit 1
fi
