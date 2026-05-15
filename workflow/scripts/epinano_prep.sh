#!/usr/bin/env bash

set -x

ulimit -n 65536 2>/dev/null || true

exec &>"${snakemake_log[0]}"

per_site="${snakemake_output[per_site]}"
sum_err="${snakemake_output[sum_err]}"

# Step 1: variant frequency extraction. Note: Epinano 1.2.0's
# Epinano_Variants.py uses `-R` (uppercase) for the reference and `-r`
# for fastq reads — easy to get wrong.
#
# Epinano 1.2.0 has an upstream bug in epinano_modules.py
# (print_last_consecutive_lines writes a 1-field 'Null' placeholder that
# slide_per_site_var then chokes on with `IndexError: ary[6]`). The crash
# happens AFTER df_proc has already written *.plus_strand.per.site.var.csv,
# and the downstream NanoRNAMod pipeline only consumes that file via
# Epinano_sumErr.py — slide_per_site_var's *.per_site.5mer.csv output is
# unused. So we tolerate a non-zero exit here as long as per_site exists.
set +e
python /opt/epinano/Epinano_Variants.py -s /opt/epinano/misc/sam2tsv.jar \
    -R "${snakemake_input[reference]}" \
    -b "${snakemake_input[sample_bam]}" \
    -t ${snakemake[threads]} -T t
rc=$?
set -e

if [ ! -f "${per_site}" ]; then
    echo "ERROR: Epinano_Variants.py did not produce ${per_site} (rc=$rc)" >&2
    echo "  Check the log above for stderr from Epinano_Variants.py." >&2
    exit 1
fi
if [ "$rc" -ne 0 ]; then
    echo "WARNING: Epinano_Variants.py exited ${rc} but ${per_site} was produced;" >&2
    echo "  this is the known slide_per_site_var IndexError in EpiNano 1.2.0." >&2
    echo "  Continuing because the downstream sumErr step only needs per_site." >&2
fi

# Step 2: per-site error summary.
python /opt/epinano/misc/Epinano_sumErr.py --kmer 0 \
    --file "${per_site}" \
    --out "${sum_err}"
