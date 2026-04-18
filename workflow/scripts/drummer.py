import os
import sys
from snakemake.shell import shell
# os.makedirs(os.path.dirname(snakemake.log[0]),exist_ok=True)
# log = snakemake.log_fmt_shell(stdout=True, stderr=True)

control_bam=snakemake.input.control_bam
native_bam=snakemake.input.native_bam
region=snakemake.input.region
reference=snakemake.params.reference
output=snakemake.output.outdir






control_bam = os.path.abspath(control_bam)
native_bam = os.path.abspath(native_bam)
region = os.path.abspath(region)
reference = os.path.abspath(reference)
output = os.path.abspath(output)

os.chdir('/opt/DRUMMER')

shell("python3 DRUMMER.py -r {reference} "
      "-l {region} "
      " -c {control_bam} "
      " -t {native_bam} "
      " -o {output} "
      " {snakemake.params.extra} "
      " -a isoform  ")

# Empty-output safety net: DRUMMER may succeed (exit 0) but produce no
# summary.txt when depth/quality filters reject all sites. Guarantee the
# output directory exists so Snakemake won't fail on the missing
# directory() output; format.py detects the absence of summary files and
# emits an empty sentinel downstream.
os.makedirs(output, exist_ok=True)

