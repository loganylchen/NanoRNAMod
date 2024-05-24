import os
import sys
from snakemake.shell import shell


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

