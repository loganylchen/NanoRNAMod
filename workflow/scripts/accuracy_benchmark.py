"""
Compute per-tool accuracy metrics against a ground-truth modification set.

Truth set TSV format (tab-separated, with header):
    transcript        - transcript identifier (must match tool output)
    position          - 0-based position on transcript
    modification_type - e.g. "m6A", "Psi", "m1Psi"

Tool result TSV format (standard NanoRNAMod output):
    transcript  - transcript identifier
    position    - 0-based position
    (other tool-specific columns)

Output columns:
    tool, modification_type, precision, recall, f1, tp, fp, fn, total_truth, total_predicted

Configuration (from snakemake.params):
    window  - positional tolerance in nucleotides (0 = exact match)
               Positions within ±window nt of a truth site count as a hit.
"""

import os
import pandas as pd

result_files = snakemake.input.results    # list of *_results.tsv paths
truth_set_path = snakemake.input.truth_set
output_file = snakemake.output[0]
window = int(snakemake.params.window)

# ── Load truth set ─────────────────────────────────────────────────────────
truth = pd.read_csv(truth_set_path, sep='\t')
required_cols = {"transcript", "position", "modification_type"}
assert required_cols.issubset(truth.columns), (
    f"Truth set must have columns: {required_cols}. Found: {set(truth.columns)}"
)

# ── Helper: extract tool name from result file path ────────────────────────
def tool_from_path(path):
    """
    Path pattern: .../modifications/{tool}/{comp_or_sample}/{tool}_results.tsv
    """
    parts = path.split(os.sep)
    for i, part in enumerate(parts):
        if part == "modifications" and i + 1 < len(parts):
            return parts[i + 1]
    return os.path.basename(os.path.dirname(os.path.dirname(path)))

# ── Load all tool results ──────────────────────────────────────────────────
tool_dfs = {}
for f in result_files:
    tool = tool_from_path(f)
    try:
        df = pd.read_csv(f, sep='\t')
        if 'transcript' in df.columns and 'position' in df.columns:
            if tool not in tool_dfs:
                tool_dfs[tool] = []
            tool_dfs[tool].append(df[['transcript', 'position']])
    except Exception as e:
        print(f"Warning: could not read {f}: {e}")

for tool in tool_dfs:
    tool_dfs[tool] = pd.concat(tool_dfs[tool], ignore_index=True).drop_duplicates()

# ─────────────────────────────────────────────────────────────────────────────
# Compute accuracy metrics for each tool and modification type
# ─────────────────────────────────────────────────────────────────────────────

records = []
mod_types = truth['modification_type'].unique()

for tool, pred_df in tool_dfs.items():
    for mod_type in mod_types:
        truth_subset = truth[truth['modification_type'] == mod_type].copy()
        if truth_subset.empty:
            continue

        # Track which predicted sites matched truth sites
        matched_pred_indices = set()

        # For each truth site, check if there's a predicted site within window
        tp = 0
        for _, truth_row in truth_subset.iterrows():
            tx = truth_row['transcript']
            pos = truth_row['position']

            # Find predictions within window on same transcript
            matches = pred_df[
                (pred_df['transcript'] == tx) &
                (pred_df['position'].between(pos - window, pos + window))
            ]

            if not matches.empty:
                tp += 1
                matched_pred_indices.update(matches.index.tolist())

        fn = len(truth_subset) - tp
        fp = len(pred_df) - len(matched_pred_indices)

        total_predicted = len(pred_df)
        total_truth = len(truth_subset)

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = (2 * precision * recall / (precision + recall)
              if (precision + recall) > 0 else 0)

        records.append({
            "tool": tool,
            "modification_type": mod_type,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "tp": tp, "fp": fp, "fn": fn,
            "total_truth": total_truth,
            "total_predicted": total_predicted,
        })

# Write output
if records:
    out_df = pd.DataFrame(records).sort_values(["modification_type", "f1"], ascending=[True, False])
else:
    out_df = pd.DataFrame(columns=["tool", "modification_type", "precision", "recall",
                                   "f1", "tp", "fp", "fn", "total_truth", "total_predicted"])
out_df.to_csv(output_file, sep='\t', index=False)
