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

# Epinano_DiffErr.R cleanup() expects a "position" column, but
# Epinano_Variants.py (v1.2+) outputs "pos". Add "position" as a copy of "pos".
add_position_col() {
    local src="$1" dst="$2"
    # pos is column 2 in: #Ref,pos,base,strand,cov,...
    awk -F',' -v OFS=',' 'NR==1{print $0,"position"} NR>1{print $0,$2}' \
        "${src}" > "${dst}"
}

control_fixed="${output_dir}/$(basename "${control}" .csv)_fixed.csv"
native_fixed="${output_dir}/$(basename "${native}" .csv)_fixed.csv"

echo "Preprocessing: adding 'position' column for Epinano_DiffErr.R compatibility"
add_position_col "${control}" "${control_fixed}"
add_position_col "${native}" "${native_fixed}"

echo "Running Epinano_DiffErr.R"
echo "  control: ${control_fixed}"
echo "  native:  ${native_fixed}"
echo "  prefix:  ${prefix}"
echo "  threads: ${threads}"
echo "  extra:   ${extra}"

Rscript /opt/epinano/Epinano_DiffErr.R \
    -t "${threads}" \
    -k "${control_fixed}" \
    -w "${native_fixed}" \
    -o "${prefix}" \
    ${extra} \
    1>"${snakemake_log[stdout]}" 2>"${snakemake_log[stderr]}"

# Clean up temp fixed files
rm -f "${control_fixed}" "${native_fixed}"

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
