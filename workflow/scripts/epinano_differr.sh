#!/usr/bin/env bash

set -e

output_dir=$(dirname "${snakemake_output[results]}")
expected_output="${snakemake_output[results]}"
control="${snakemake_input[control]}"
native="${snakemake_input[native]}"
prefix="${snakemake_params[prefix]}"
extra="${snakemake_params[extra]}"

mkdir -p "${output_dir}"

# Epinano_DiffErr.R (v1.2.0) arguments:
#   -k: knockout/unmodified sample (control)
#   -w: wildtype/modified sample (native)
#   -o: output prefix
#   -t: z-score threshold (NOT threads)
#   -c: minimum coverage
#   -f: feature column to analyze (e.g. sum_err)
#   -d: minimum deviance
#   -p: generate plots (optional)
# All of -t, -c, -f, -d are passed via ${extra} from config.

echo "Running Epinano_DiffErr.R"
echo "  control: ${control}"
echo "  native:  ${native}"
echo "  prefix:  ${prefix}"
echo "  extra:   ${extra}"

Rscript /opt/Epinano/Epinano_DiffErr.R \
    -k "${control}" \
    -w "${native}" \
    -o "${prefix}" \
    ${extra} \
    1>"${snakemake_log[0]}" 2>&1

echo "Files produced in ${output_dir}:"
ls -lh "${output_dir}/" 2>/dev/null || echo "  (directory empty)"

# If the expected file exists, we're done
if [ -f "${expected_output}" ]; then
    echo "Output found: ${expected_output}"
    exit 0
fi

# EpiNano may name the file differently — find any *.prediction.csv in the dir
found=$(find "${output_dir}" -maxdepth 1 -name "*.prediction.csv" -print -quit 2>/dev/null)
if [ -n "${found}" ]; then
    echo "Renaming ${found} -> ${expected_output}"
    mv "${found}" "${expected_output}"
else
    echo "ERROR: Epinano_DiffErr.R produced no *.prediction.csv in ${output_dir}"
    echo "Check log: ${snakemake_log[0]}"
    exit 1
fi
