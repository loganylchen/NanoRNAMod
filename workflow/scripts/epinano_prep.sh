#!/usr/bin/env bash

set -x
set -e

exec &>"${snakemake_log[0]}"

per_site="${snakemake_output[per_site]}"
sum_err="${snakemake_output[sum_err]}"

# Step 1: variant frequency extraction. Note: Epinano 1.2.0's
# Epinano_Variants.py uses `-R` (uppercase) for the reference and `-r`
# for fastq reads — easy to get wrong.
python /opt/epinano/Epinano_Variants.py -s /opt/epinano/misc/sam2tsv.jar \
    -R "${snakemake_input[reference]}" \
    -b "${snakemake_input[sample_bam]}" \
    -t ${snakemake[threads]} -T t

# Epinano_Variants.py sometimes exits 0 without producing the expected
# .plus_strand.per.site.var.csv (e.g. argparse validation failures call
# bare exit(), or df_proc skips the write when df_is_not_empty is false).
# Fail loudly here so the real error surfaces instead of masking it.
if [ ! -f "${per_site}" ]; then
    echo "ERROR: Epinano_Variants.py exited 0 but did not produce ${per_site}" >&2
    echo "  Check the log above for stderr from Epinano_Variants.py." >&2
    exit 1
fi

# Step 2: per-site error summary.
python /opt/epinano/misc/Epinano_sumErr.py --kmer 0 \
    --file "${per_site}" \
    --out "${sum_err}"
