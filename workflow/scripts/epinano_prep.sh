#!/usr/bin/env bash

#minimap2 ${snakemake_params[opts]} -t ${snakemake[threads]} "${snakemake_input[reference]}" \
#    "${snakemake_input[0]}" > "${snakemake_output[0]}" 2>> "${snakemake_log[0]}"
set -x
set -e

python /opt/Epinano/Epinano_Variants.py  -s /opt/Epinano/misc/sam2tsv.jar  -R "${snakemake_input[reference]}" \
-b "${snakemake_input[sample_bam]}" -n {snakemake[threads]} -T t 2>"${snakemake_log[0]}"

python /opt/Epinano/misc/Epinano_sumErr.py --kmer 0 --file "${snakemake_output[per_site]}" \
    --out "${snakemake_output[sum_err]}" 2>> "${snakemake_log[0]}"

#python /opt/Epinano/misc/Slide_Variants.py "${snakemake_output[per_site]}" 5
#python /opt/Epinano/Epinano_Predict.py  \
#    --model /opt/Epinano/models/rrach.q3.mis3.del3.linear.dump \
#    --predict "${snakemake_output[kmer_5_site]}" \
#    --columns 8,13,23 \
#    --out_prefix "${snakemake_params[prefix]}"