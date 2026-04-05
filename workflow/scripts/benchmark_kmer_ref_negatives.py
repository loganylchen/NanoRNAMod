"""
Reference-based 5-mer matched negative site generation.

Scans the reference transcriptome FASTA for all positions whose 5-mer context
matches a 5-mer found at a positive (ground-truth) site. Positions NOT in the
positive set become negatives. This controls for sequence-context bias when the
ground truth contains only positives.

Input:
    truth_set: TSV with columns transcript, position, modification_type
    transcriptome_fasta: reference FASTA (transcript sequences)

Output:
    negatives: TSV with columns transcript, position, kmer, label ("-")

Params:
    window: positions within this distance of a positive are excluded
"""

import pandas as pd

truth_path = snakemake.input.truth_set
fasta_path = snakemake.input.transcriptome_fasta
output_path = snakemake.output.negatives
window = int(snakemake.params.window)


def read_fasta(path):
    """Parse a FASTA file into a dict of {name: sequence}."""
    import gzip

    sequences = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        name = None
        parts = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            sequences[name] = "".join(parts)
    return sequences


# Load truth set
truth_df = pd.read_csv(truth_path, sep="\t")
if "label" in truth_df.columns:
    truth_pos = truth_df[truth_df["label"] != "-"].copy()
else:
    truth_pos = truth_df.copy()

truth_pos["position"] = pd.to_numeric(truth_pos["position"], errors="coerce")
truth_pos = truth_pos.dropna(subset=["position"])
truth_pos["position"] = truth_pos["position"].astype(int)

print(f"Loaded {len(truth_pos)} positive sites from truth set")

# Build set of positive (transcript, position) for exclusion
positive_set = set(zip(truth_pos["transcript"], truth_pos["position"]))

# Build dict of positive positions per transcript for window exclusion
positive_by_tx = {}
for tx, pos in positive_set:
    positive_by_tx.setdefault(tx, []).append(pos)

# Load reference transcriptome
fasta = read_fasta(fasta_path)
print(f"Loaded {len(fasta)} transcripts from reference FASTA")

# Extract 5-mers at positive sites from reference
positive_kmers = set()
flanking = 2  # 5-mer: center +-2

for _, row in truth_pos.iterrows():
    tx = row["transcript"]
    pos = int(row["position"])
    if tx not in fasta:
        continue
    seq = fasta[tx]
    start = pos - flanking
    end = pos + flanking + 1
    if start < 0 or end > len(seq):
        continue
    kmer = seq[start:end].upper()
    if len(kmer) == 5 and "N" not in kmer:
        positive_kmers.add(kmer)

print(f"Extracted {len(positive_kmers)} unique 5-mers from positive sites")

if not positive_kmers:
    print("WARNING: No 5-mers extracted. Writing empty output.")
    pd.DataFrame(columns=["transcript", "position", "kmer", "label"]).to_csv(
        output_path, sep="\t", index=False
    )
    raise SystemExit(0)

# Scan all transcripts for matching 5-mers
negatives = []
for tx_name, seq in fasta.items():
    seq = seq.upper()
    tx_positives = positive_by_tx.get(tx_name, [])

    for i in range(flanking, len(seq) - flanking):
        kmer = seq[i - flanking : i + flanking + 1]
        if len(kmer) != 5:
            continue
        if kmer not in positive_kmers:
            continue

        # Exclude if within window of any positive site on this transcript
        too_close = False
        for pp in tx_positives:
            if abs(i - pp) <= window:
                too_close = True
                break
        if too_close:
            continue

        negatives.append((tx_name, i, kmer))

neg_df = pd.DataFrame(negatives, columns=["transcript", "position", "kmer"])
neg_df["label"] = "-"
neg_df = neg_df.drop_duplicates(subset=["transcript", "position"])

print(f"Generated {len(neg_df)} reference-based 5-mer matched negative sites")
print(f"  Across {neg_df['transcript'].nunique()} transcripts")

neg_df.to_csv(output_path, sep="\t", index=False)
