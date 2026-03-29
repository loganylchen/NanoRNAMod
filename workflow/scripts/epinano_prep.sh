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

# Epinano_Variants.py (v1.2.0) derives output prefix from the BAM path:
#   prefix = bam_file.replace('.bam','')
#   output = prefix + '.plus_strand.per.site.var.csv'
# So the output lands in the same directory as the BAM file.
# We cd to the target directory as a safety measure for any version differences.
( cd "${directory}" && \
  python /opt/epinano/Epinano_Variants.py -R "${ref_abs}" \
      -s /opt/epinano/misc/sam2tsv.jar \
      -b "${bam_abs}" -t ${snakemake[threads]} ${extra} ) \
    2>&1 | tee -a "${snakemake_log[0]}"

# Check pipeline exit status (not tee's exit status)
epinano_exit=${PIPESTATUS[0]}
if [ ${epinano_exit} -ne 0 ]; then
    echo "ERROR: Epinano_Variants.py failed with exit code ${epinano_exit}" | tee -a "${snakemake_log[0]}"
    exit ${epinano_exit}
fi

# Verify output exists
if [ ! -f "${expected_output}" ]; then
    echo "ERROR: EpiNano did not produce expected output: ${expected_output}" | tee -a "${snakemake_log[0]}"
    echo "Files in ${directory}:" | tee -a "${snakemake_log[0]}"
    ls -la "${directory}"/*.per.site*.csv 2>/dev/null | tee -a "${snakemake_log[0]}" || echo "  (no .per.site*.csv files found)" | tee -a "${snakemake_log[0]}"
    exit 1
fi
