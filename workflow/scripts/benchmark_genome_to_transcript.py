"""
Convert genome-coordinate truth set to transcript coordinates using GTF exon annotations.

For each genome position in the truth set, finds all transcripts whose exons
overlap that position and computes the corresponding transcript coordinate.
One genome position may map to multiple transcripts (all rows emitted).
Positions falling in introns are skipped.

Input:
    truth_set: TSV with genome-level columns (chrom + position) plus modification_type, label
    gtf: GTF file with exon features for coordinate mapping

Output:
    converted: TSV with transcript-level columns (transcript, position) plus
               original chrom, genome_position for traceability
"""

import csv
import sys
from collections import defaultdict

truth_path = snakemake.input.truth_set
gtf_path = snakemake.input.gtf
output_path = snakemake.output.converted

log_path = snakemake.log[0] if snakemake.log else None
log_fh = open(log_path, "w") if log_path else sys.stderr


def log(msg):
    print(msg, file=log_fh)


# --- Step 1: Parse GTF for exon coordinates per transcript ---

# transcript_id -> { "chrom": str, "strand": str, "exons": [(start, end), ...] }
transcripts = defaultdict(lambda: {"chrom": None, "strand": None, "exons": []})

log("Parsing GTF...")
with open(gtf_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        if fields[2] != "exon":
            continue

        chrom = fields[0]
        start = int(fields[3]) - 1  # GTF is 1-based, convert to 0-based
        end = int(fields[4])  # GTF end is inclusive, but as 0-based half-open this is correct
        strand = fields[6]

        # Parse transcript_id from attributes
        attrs = fields[8]
        tid = None
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("transcript_id"):
                # Handle both: transcript_id "X" and transcript_id X
                tid = attr.split('"')[1] if '"' in attr else attr.split()[1]
                break

        if tid is None:
            continue

        rec = transcripts[tid]
        rec["chrom"] = chrom
        rec["strand"] = strand
        rec["exons"].append((start, end))

# Sort exons by genomic position for each transcript
for tid, rec in transcripts.items():
    rec["exons"].sort(key=lambda x: x[0])

log(f"Parsed {len(transcripts)} transcripts from GTF")

# Build chrom -> list of (start, end, tid) index for fast lookup
chrom_index = defaultdict(list)
for tid, rec in transcripts.items():
    for start, end in rec["exons"]:
        chrom_index[rec["chrom"]].append((start, end, tid))

# Sort each chrom's exon list by start position
for chrom in chrom_index:
    chrom_index[chrom].sort(key=lambda x: x[0])


def genome_to_transcript_pos(tid, genome_pos, rec):
    """Convert a 0-based genome position to a 0-based transcript position.

    Returns None if the position falls in an intron (not covered by any exon).
    For minus-strand transcripts, the position is measured from the 3' end
    (transcript start) which corresponds to the last exon in genome order.
    """
    exons = rec["exons"]
    strand = rec["strand"]

    # Calculate transcript length and cumulative offset
    transcript_offset = 0
    hit_offset = None

    for ex_start, ex_end in exons:
        if ex_start <= genome_pos < ex_end:
            hit_offset = transcript_offset + (genome_pos - ex_start)
        transcript_offset += ex_end - ex_start

    if hit_offset is None:
        return None

    transcript_length = transcript_offset

    if strand == "-":
        return transcript_length - 1 - hit_offset
    return hit_offset


# --- Step 2: Load truth set and detect columns ---

import pandas as pd

truth_df = pd.read_csv(truth_path, sep="\t")

# Normalize chrom/position columns
chrom_candidates = ["chrom", "chr", "chromosome", "seqname"]
pos_candidates = ["position", "pos", "start", "genomic_position", "genome_position"]

chrom_col = None
for c in chrom_candidates:
    if c in truth_df.columns:
        chrom_col = c
        break

pos_col = None
for c in pos_candidates:
    if c in truth_df.columns:
        pos_col = c
        break

if chrom_col is None or pos_col is None:
    raise ValueError(
        f"Truth set must have chrom column (tried {chrom_candidates}) "
        f"and position column (tried {pos_candidates}). "
        f"Found columns: {list(truth_df.columns)}"
    )

log(f"Using chrom column: '{chrom_col}', position column: '{pos_col}'")
log(f"Truth set has {len(truth_df)} rows")

# --- Step 3: Convert each position ---

rows_out = []
skipped_no_chrom = 0
skipped_intron = 0
converted = 0

for _, row in truth_df.iterrows():
    chrom = str(row[chrom_col])
    genome_pos = int(row[pos_col])

    if chrom not in chrom_index:
        skipped_no_chrom += 1
        continue

    # Find all transcripts with exons overlapping this position
    seen_tids = set()
    for ex_start, ex_end, tid in chrom_index[chrom]:
        if ex_end <= genome_pos:
            continue
        if ex_start > genome_pos:
            break
        if tid in seen_tids:
            continue
        seen_tids.add(tid)

        rec = transcripts[tid]
        tx_pos = genome_to_transcript_pos(tid, genome_pos, rec)
        if tx_pos is None:
            skipped_intron += 1
            continue

        out_row = {
            "transcript": tid,
            "position": tx_pos,
            "chrom": chrom,
            "genome_position": genome_pos,
        }
        # Carry over additional columns
        for col in truth_df.columns:
            if col not in (chrom_col, pos_col) and col not in out_row:
                out_row[col] = row[col]

        rows_out.append(out_row)
        converted += 1

    if not seen_tids:
        skipped_intron += 1

log(f"Converted {converted} position-transcript pairs")
log(f"Skipped {skipped_no_chrom} rows (chrom not in GTF)")
log(f"Skipped {skipped_intron} rows (intron / no overlapping transcript)")

# --- Step 4: Write output ---

out_df = pd.DataFrame(rows_out)

if out_df.empty:
    log("WARNING: No positions could be converted. Output will be empty.")
    # Write header-only file so downstream rules don't fail on missing file
    out_df = pd.DataFrame(
        columns=["transcript", "position", "chrom", "genome_position", "modification_type", "label"]
    )

# Ensure standard column order: transcript, position first
priority_cols = ["transcript", "position", "modification_type", "label", "chrom", "genome_position"]
other_cols = [c for c in out_df.columns if c not in priority_cols]
col_order = [c for c in priority_cols if c in out_df.columns] + other_cols
out_df = out_df[col_order]

out_df.to_csv(output_path, sep="\t", index=False)
log(f"Wrote {len(out_df)} rows to {output_path}")

if log_fh is not sys.stderr:
    log_fh.close()
