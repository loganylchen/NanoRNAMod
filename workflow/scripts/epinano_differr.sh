#!/usr/bin/env bash

set -e

output_dir=$(dirname "${snakemake_output[results]}")
expected_output="${snakemake_output[results]}"
control="${snakemake_input[control]}"
native="${snakemake_input[native]}"
prefix="${snakemake_params[prefix]}"
threads="${snakemake[threads]}"

# Strip any -t/--threads flags from extra params to avoid duplicates
extra=$(echo "${snakemake_params[extra]}" | sed -E 's/-t[[:space:]]+[0-9]+//g')

mkdir -p "${output_dir}"

echo "Running Epinano_DiffErr.R"
echo "  control: ${control}"
echo "  native:  ${native}"
echo "  prefix:  ${prefix}"
echo "  threads: ${threads}"
echo "  extra:   ${extra}"

Rscript /opt/epinano/Epinano_DiffErr.R \
    -t "${threads}" \
    -k "${control}" \
    -w "${native}" \
    -o "${prefix}" \
    ${extra} \
    1>"${snakemake_log[stdout]}" 2>"${snakemake_log[stderr]}"

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
    echo "Check logs: ${snakemake_log[stderr]}"
    exit 1
fi
