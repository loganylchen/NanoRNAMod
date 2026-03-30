#!/usr/bin/env bash

set -x
set -e

python /opt/epinano/Epinano_Variants.py -s /opt/epinano/misc/sam2tsv.jar \
    -R "${snakemake_input[reference]}" \
    -b "${snakemake_input[sample_bam]}" \
    -t ${snakemake[threads]} -T t \
    2>"${snakemake_log[0]}"

python /opt/epinano/misc/Epinano_sumErr.py --kmer 0 \
    --file "${snakemake_output[per_site]}" \
    --out "${snakemake_output[sum_err]}" \
    2>>"${snakemake_log[0]}"
