#!/usr/bin/env bash

#minimap2 ${snakemake_params[opts]} -t ${snakemake[threads]} "${snakemake_input[reference]}" \
#    "${snakemake_input[0]}" > "${snakemake_output[0]}" 2>> "${snakemake_log[0]}"
set -x
set -e

python /opt/epinano/Epinano_Variants.py   -r "${snakemake_input[reference]}" \
-b "${snakemake_input[sample_bam]}" -c ${snakemake[threads]} -o "${snakemake_output[per_site]}" 2>"${snakemake_log[0]}"

ls -lh "$(dirname "${snakemake_output[per_site]}")"


# python /opt/epinano/misc/Epinano_sumErr.py --kmer 0 --file "${snakemake_output[per_site]}" \
#     --out "${snakemake_output[sum_err]}" 2>> "${snakemake_log[0]}"

#python /opt/epinano/misc/Slide_Variants.py "${snakemake_output[per_site]}" 5
#python /opt/epinano/Epinano_Predict.py  \
#    --model /opt/epinano/models/rrach.q3.mis3.del3.linear.dump \
#    --predict "${snakemake_output[kmer_5_site]}" \
#    --columns 8,13,23 \
#    --out_prefix "${snakemake_params[prefix]}"