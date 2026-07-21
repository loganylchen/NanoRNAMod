import os
import sys
import glob
import shutil
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

# DRUMMER (isoform mode) creates one set of BAM/analysis files *per transcript*
# inside its -o directory. For a full transcriptome that is hundreds of
# thousands of tiny files. Writing them directly onto a shared/network
# filesystem makes DRUMMER severely I/O-bound (observed ~0.2 effective cores,
# multi-day runtimes) and burns an enormous number of inodes.
#
# Instead we run DRUMMER against node-local scratch and copy back only the
# per-comparison summary.txt files -- the sole output consumed downstream by
# scripts/format.py::format_drummer (which globs "<control>*-<native>*/summary.txt").
# Scratch defaults to the Snakemake `tmpdir` resource, then $TMPDIR, then /tmp.
scratch_root = getattr(snakemake.resources, "tmpdir", None) \
    or os.environ.get("TMPDIR") or "/tmp"
local_out = os.path.join(scratch_root, "drummer_" + os.path.basename(output.rstrip("/")))
shutil.rmtree(local_out, ignore_errors=True)
os.makedirs(local_out, exist_ok=True)

os.chdir('/opt/DRUMMER')

try:
    shell("python3 DRUMMER.py -r {reference} "
          "-l {region} "
          " -c {control_bam} "
          " -t {native_bam} "
          " -o {local_out} "
          " {snakemake.params.extra} "
          " -a isoform  ")

    # Copy back only the per-comparison summary.txt files (one per
    # "<control>-<native>" sub-directory). The per-transcript intermediates
    # stay on node-local scratch and are discarded below.
    os.makedirs(output, exist_ok=True)
    for summary in glob.glob(os.path.join(local_out, "*", "summary.txt")):
        dst_dir = os.path.join(output, os.path.basename(os.path.dirname(summary)))
        os.makedirs(dst_dir, exist_ok=True)
        shutil.copy2(summary, dst_dir)
finally:
    shutil.rmtree(local_out, ignore_errors=True)

# Empty-output safety net: DRUMMER may succeed (exit 0) but produce no
# summary.txt when depth/quality filters reject all sites. Guarantee the
# output directory exists so Snakemake won't fail on the missing
# directory() output; format.py detects the absence of summary files and
# emits an empty sentinel downstream.
os.makedirs(output, exist_ok=True)
