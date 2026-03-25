"""
Compute per-tool accuracy metrics against a ground-truth modification set.

Truth set TSV format (tab-separated, with header):
    transcript        - transcript identifier (must match tool output)
    position          - 0-based position on transcript
    modification_type - e.g. "m6A", "Psi", "m1Psi"
    label             - (optional) "-" for negative sites, any other value for positive
                       If absent, all sites are treated as positive.
                       Negative sites will be auto-inferred if not explicitly provided.

Tool result TSV format (standard NanoRNAMod output):
    transcript  - transcript identifier (or id, chrom, ref_id depending on tool)
    position    - 0-based position (or pos, start, start_loc depending on tool)
    (other tool-specific columns including optional scores)

Output columns:
    tool, modification_type, window, precision, recall, f1, tp, fp, fn, tn,
    specificity, mcc, auprc, auroc, called_sites,
    total_truth, total_predicted, total_negative

Configuration (from snakemake.params):
    window  - positional tolerance in nucleotides (0 = exact match)
               Can be a single integer or a list of integers for multi-window evaluation.
               Positions within ±window nt of a truth site count as a hit.
"""

import os
import ast
import pandas as pd
import numpy as np

result_files = snakemake.input.results
truth_set_path = snakemake.input.truth_set

# Named outputs from Snakemake rule
aggregated_files = {
    "by_mod_type": snakemake.output.summary,           # accuracy_summary.tsv
    "overall": snakemake.output.overall,               # accuracy_summary_overall.tsv
    "by_comparison": snakemake.output.by_comparison,   # accuracy_summary_by_comparison.tsv
    "by_negative_type": snakemake.output.by_negative_type,  # accuracy_summary_by_negative_type.tsv
}

# Derive benchmark directory from first output file
benchmark_dir = os.path.dirname(snakemake.output.summary)

window_param = snakemake.params.window

# Support both single int and list of windows
if isinstance(window_param, int):
    windows = [window_param]
elif isinstance(window_param, list):
    windows = window_param
else:
    try:
        parsed = ast.literal_eval(str(window_param))
        if isinstance(parsed, list):
            windows = parsed
        elif isinstance(parsed, int):
            windows = [parsed]
        else:
            windows = [int(window_param)]
    except (ValueError, SyntaxError):
        windows = [int(window_param)]

# ── Load truth set ─────────────────────────────────────────────────────────
truth = pd.read_csv(truth_set_path, sep='\t')
required_cols = {"transcript", "position", "modification_type"}
assert required_cols.issubset(truth.columns), (
    f"Truth set must have columns: {required_cols}. Found: {set(truth.columns)}"
)

# Check for label column (optional)
has_label = 'label' in truth.columns

if has_label:
    # Split truth into positive and negative sites
    truth_pos = truth[truth['label'] != '-'].copy()
    truth_neg = truth[truth['label'] == '-'].copy()
else:
    # All sites are positive; negatives will be auto-inferred
    truth_pos = truth.copy()
    truth_neg = pd.DataFrame(columns=truth.columns)

# ── Helper: extract tool name from result file path ────────────────────────
def tool_from_path(path):
    """Extract tool name from path like .../modifications/{tool}/{comp}/{tool}_results.tsv
    or .../results/{tool}/{comp}/transcript_mod_results.csv"""
    parts = path.split(os.sep)
    for i, part in enumerate(parts):
        if part == "modifications" and i + 1 < len(parts):
            return parts[i + 1]
        if part == "results" and i + 1 < len(parts):
            # Check if next part looks like a tool name (not 'alignments', 'fastq', etc.)
            potential_tool = parts[i + 1]
            non_tool_dirs = ['alignments', 'fastq', 'blow5', 'modifications', 'benchmarks', 'viz']
            if potential_tool not in non_tool_dirs:
                return potential_tool
    return os.path.basename(os.path.dirname(os.path.dirname(path)))


def comparison_from_path(path):
    """Extract comparison name from path like .../modifications/{tool}/{comp}/{tool}_results.tsv
    or .../results/{tool}/{comp}/transcript_mod_results.csv"""
    parts = path.split(os.sep)
    for i, part in enumerate(parts):
        if part == "modifications" and i + 2 < len(parts):
            return parts[i + 2]
        if part == "results" and i + 2 < len(parts):
            # Check if this looks like a comparison directory
            potential_comp = parts[i + 2]
            if '_' in potential_comp:  # Comparison names typically have underscore (e.g., native_0_control_0)
                return potential_comp
    return "unknown"

# ── Helper: normalize column names for transcript and position ────────────
def normalize_columns(df):
    """Map various column names to standard 'transcript' and 'position'."""
    col_mapping = {}

    # Transcript column variants
    transcript_cols = ['transcript', 'id', 'ref_id', 'chrom']
    for col in transcript_cols:
        if col in df.columns:
            col_mapping[col] = 'transcript'
            break

    # Position column variants
    position_cols = ['position', 'pos', 'start', 'start_loc', 'transcript_pos', 'transcript_loc', 'genomic_pos']
    for col in position_cols:
        if col in df.columns and col not in col_mapping.values():
            col_mapping[col] = 'position'
            break

    if col_mapping:
        df = df.rename(columns=col_mapping)

    return df

# ── Helper: detect score column in dataframe ───────────────────────────────
def detect_score_column(df, tool_name=None):
    """Detect the most likely score column for ranking predictions.

    Score column priority per tool:
    - P-value columns: lower is better (converted to -log10 for ranking)
    - Probability/score columns: higher is better

    Column names are matched using:
    1. Exact match
    2. Case-insensitive match
    3. Substring match (for columns with prefixes/suffixes like pval_CASE_vs_CONTROL)
    """
    # Tool-specific score column patterns (in order of preference)
    # Based on documented tool outputs from tools_comparison.md and actual output files
    # Priority: p-value/adjusted p-value first (for comparative tools),
    #           then probability/score (for ML tools)
    tool_score_patterns = {
        # Signal-based comparative tools
        # xPore raw output: id,position,kmer,diff_mod_rate_CASE_vs_CONTROL,pval_CASE_vs_CONTROL,z_score_CASE_vs_CONTROL,...
        # format.py preserves original column names
        'xpore': ['pval_CASE_vs_CONTROL', 'z_score_CASE_vs_CONTROL', 'diff_mod_rate_CASE_vs_CONTROL',
                  'p_value', 'adjusted_p_value', 'statistic', 'pval', 'padj'],
        # Nanocompore raw output: pos,chr,genomicPos,ref_id,strand,ref_kmer,GMM_logit_pvalue,KS_dwell_pvalue,KS_intensity_pvalue,...
        'nanocompore': ['GMM_logit_pvalue', 'KS_dwell_pvalue', 'KS_intensity_pvalue', 'Logit_LOR',
                        'p_value_glm', 'p_value_ks', 'p_value'],
        # Baleen raw output: transcript,position,kmer,homopolymer,...,p_value,padj,...
        # Full columns: transcript,position,kmer,homopolymer,control_depth,native_depth,control_LOD,native_LOD,
        #               control_normalized_LOD,native_normalized_LOD,nLOD_diff,control_mod_ratio,native_mod_ratio,
        #               mod_ratio_diff,mod_ratio_fc,z_test,p_value,padj
        'baleen': ['p_value', 'padj', 'z_test', 'mod_ratio_diff', 'mod_ratio_fc', 'nLOD_diff',
                   'native_mod_ratio', 'native_normalized_LOD', 'native_LOD',
                   'pvalue', 'effect_size', 'score'],
        # pyBaleen raw output: contig,position,kmer,mod_ratio,ci_low,ci_high,pvalue,padj,effect_size,...
        'pybaleen': ['pvalue', 'padj', 'effect_size', 'p_value', 'adjusted_p_value', 'stoichiometry'],

        # Error-based comparative tools
        # DiffErr raw output (BED, after format.py): -log10 P value, -log10 FDR, odds ratio, G statistic
        'differr': ['-log10 P value', '-log10 FDR', 'odds ratio', 'G statistic', 'score'],
        # DRUMMER raw output: transcript_id,reference_base,transcript_pos,...,OR_padj,...,G_test,G_padj,...
        'drummer': ['OR_padj', 'G_padj', 'G_test', 'log2_(OR)', 'p_value', 'pvalue', 'padj', 'odds_ratio'],
        # ELIGOS2 raw output: chrom,start_loc,end_loc,strand,name,ref,...,pval,adjPval,...
        'eligos2': ['pval', 'adjPval', 'p_value', 'q_value', 'error_rate', 'pvalue', 'qvalue'],
        # EpiNano raw output (after format.py): chrom,pos,ref,strand,...,delta_sum_err,z_scores,z_score_prediction
        'epinano': ['z_scores', 'delta_sum_err', 'z_score_prediction', 'score'],
        # psipore output: transcript,position,psi_score,p_value,coverage
        'psipore': ['p_value', 'psi_score', 'pvalue', 'score'],

        # Single-sample ML/DL tools
        # TandemMod output: transcript,position,modification,probability,stoichiometry
        'tandemmod': ['probability', 'stoichiometry', 'prob', 'score'],
        # DirectRM output: transcript,position,modification,probability
        'directrm': ['probability', 'prob', 'score'],
        # m6ATM output: transcript,position,probability,stoichiometry
        'm6atm': ['probability', 'stoichiometry', 'prob', 'score'],
        # Rnano output: transcript,position,modification,probability
        'rnano': ['probability', 'prob', 'score'],
        # NanoPSU output: transcript,position,psi_probability,coverage
        'nanopsu': ['psi_probability', 'probability', 'pvalue', 'score'],
        # NanoMUD output: transcript,position,modification,probability,coverage
        'nanomud': ['probability', 'pvalue', 'p_value', 'score'],
        # Penguin output: transcript,position,psi_probability,score
        'penguin': ['psi_probability', 'probability', 'score', 'pvalue'],
    }

    def match_column(df_columns, pattern):
        """Match a pattern against dataframe columns using multiple strategies."""
        pattern_lower = pattern.lower().replace(' ', '_')

        # 1. Exact match
        if pattern in df_columns:
            return pattern

        # 2. Case-insensitive exact match
        for col in df_columns:
            if col.lower().replace(' ', '_') == pattern_lower:
                return col

        # 3. Substring match (pattern is contained in column name)
        for col in df_columns:
            col_lower = col.lower().replace(' ', '_')
            if pattern_lower in col_lower:
                return col

        # 4. Reverse substring match (column name is contained in pattern)
        for col in df_columns:
            col_lower = col.lower().replace(' ', '_')
            if col_lower in pattern_lower:
                return col

        return None

    df_columns = list(df.columns)

    # If tool name is provided, try tool-specific columns first
    if tool_name and tool_name in tool_score_patterns:
        for pattern in tool_score_patterns[tool_name]:
            matched = match_column(df_columns, pattern)
            if matched:
                return matched

    # Generic score columns (fallback)
    generic_patterns = [
        'p_value', 'pvalue', 'pval', '-log10_pvalue', '-log10 P value',
        'score', 'probability', 'prob', 'mod_score', 'mod_prob',
        'logit_pvalue', 'logit', 'z_score', 'z_scores'
    ]

    for pattern in generic_patterns:
        matched = match_column(df_columns, pattern)
        if matched:
            return matched

    return None

# ── Helper: compute AUPRC and AUROC ───────────────────────────────────────
def compute_ranking_metrics(pred_df, truth_pos_subset, truth_neg_subset, score_col, window=0, tool_name="unknown"):
    """
    Compute AUPRC and AUROC given predictions with scores.

    Uses binary labels based on match with ground truth sites within window.
    Score is converted to -log10(p-value) if it's a p-value column.
    """
    if score_col is None or pred_df.empty:
        return np.nan, np.nan

    # Check if score column exists and has valid values
    if score_col not in pred_df.columns:
        return np.nan, np.nan

    # Determine if this is a p-value column (needs -log10 transformation)
    # Note: columns already containing '-log10' are pre-transformed (higher = better)
    score_col_lower = score_col.lower().replace(' ', '_')
    is_log_transformed = '-log10' in score_col_lower or 'log10_pvalue' in score_col_lower
    is_pvalue_col = not is_log_transformed and any(x in score_col_lower for x in ['p_value', 'pvalue', 'pval', 'padj', 'q_value', 'qval', 'fdr', 'adj'])

    # Create labels and scores arrays
    labels = []
    scores = []
    nan_count = 0
    non_numeric_count = 0

    for _, pred_row in pred_df.iterrows():
        tx = pred_row['transcript']
        pos = pred_row['position']
        raw_score = pred_row[score_col]

        # Skip if score is NaN or non-numeric
        if pd.isna(raw_score):
            nan_count += 1
            continue
        try:
            raw_score = float(raw_score)
        except (ValueError, TypeError):
            non_numeric_count += 1
            continue

        # Convert p-value to -log10(p-value) for ranking (higher = better)
        if is_pvalue_col:
            # Clamp p-value to avoid log(0)
            pval = max(raw_score, 1e-300)
            score = -np.log10(pval)
        else:
            score = raw_score

        # Check if this prediction matches a positive truth site (within window)
        is_positive = False
        for _, truth_row in truth_pos_subset.iterrows():
            if (truth_row['transcript'] == tx and
                abs(truth_row['position'] - pos) <= window):
                is_positive = True
                break

        labels.append(1 if is_positive else 0)
        scores.append(score)

    n_valid = len(scores)
    n_positive = sum(labels)
    n_negative = n_valid - n_positive

    # Need at least one positive and one negative sample
    if len(set(labels)) < 2 or len(scores) == 0:
        return np.nan, np.nan

    try:
        from sklearn.metrics import roc_auc_score, average_precision_score
        auroc = roc_auc_score(labels, scores)
        auprc = average_precision_score(labels, scores)
        return auprc, auroc
    except Exception:
        return np.nan, np.nan

# ── Load all tool results (preserving all columns for scores) ────────────────
# Now also track comparison information for per-comparison analysis
tool_dfs = {}
tool_score_cols = {}
tool_comparisons = {}  # tool -> list of (comparison, dataframe) tuples

for f in result_files:
    tool = tool_from_path(f)
    comparison = comparison_from_path(f)

    try:
        # Handle both TSV and CSV files (baleen outputs CSV)
        if f.endswith('.csv'):
            df = pd.read_csv(f)
        else:
            df = pd.read_csv(f, sep='\t')
        df = normalize_columns(df)

        if 'transcript' in df.columns and 'position' in df.columns:
            # Detect score column for this tool (pass tool name for better detection)
            score_col = detect_score_column(df, tool)
            if score_col and tool not in tool_score_cols:
                tool_score_cols[tool] = score_col
                print(f"Detected score column for {tool}: {score_col}")

            if tool not in tool_dfs:
                tool_dfs[tool] = []
                tool_comparisons[tool] = []

            # Add comparison column to track source
            df['comparison'] = comparison
            tool_dfs[tool].append(df)
            tool_comparisons[tool].append((comparison, df))
    except Exception as e:
        print(f"Warning: could not read {f}: {e}")

for tool in tool_dfs:
    tool_dfs[tool] = pd.concat(tool_dfs[tool], ignore_index=True).drop_duplicates(
        subset=['transcript', 'position', 'comparison']
    )
    print(f"Loaded {len(tool_dfs[tool])} predictions for {tool} from {len(tool_comparisons[tool])} comparisons")

# ─────────────────────────────────────────────────────────────────────────────
# Compute accuracy metrics for each tool, modification type, and window
# With FAIR COMPARISON: sites covered by any tool get score=0 for tools that didn't call them
# ─────────────────────────────────────────────────────────────────────────────

def compute_mcc(tp, fp, fn, tn):
    """Compute Matthews Correlation Coefficient."""
    numerator = (tp * tn) - (fp * fn)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator == 0:
        return 0.0
    return numerator / denominator


def compute_ranking_metrics_with_fair_comparison(pred_df, all_covered_sites,
                                                   truth_pos_subset, truth_neg_subset,
                                                   score_col, window=0, tool_name="unknown"):
    """
    Compute AUPRC and AUROC with fair comparison.

    Fair comparison: if a site is covered by any tool but not called by this tool,
    treat it as a negative prediction with score=0.

    Parameters:
        pred_df: DataFrame with this tool's predictions
        all_covered_sites: Set of (transcript, position) tuples covered by ANY tool
        truth_pos_subset: DataFrame of positive truth sites
        truth_neg_subset: DataFrame of negative truth sites (if available)
        score_col: Column name for score
        window: Position tolerance window
        tool_name: Tool name for debugging
    """
    if pred_df.empty and not all_covered_sites:
        return np.nan, np.nan

    # Determine if score column needs transformation
    is_pvalue_col = False
    is_log_transformed = False
    if score_col is not None and score_col in pred_df.columns:
        score_col_lower = score_col.lower().replace(' ', '_')
        is_log_transformed = '-log10' in score_col_lower or 'log10_pvalue' in score_col_lower
        is_pvalue_col = not is_log_transformed and any(
            x in score_col_lower for x in ['p_value', 'pvalue', 'pval', 'padj', 'q_value', 'qval', 'fdr', 'adj']
        )

    # Build prediction lookup for this tool
    pred_lookup = {}  # (tx, pos) -> score
    if score_col is not None and score_col in pred_df.columns:
        for _, row in pred_df.iterrows():
            tx = row['transcript']
            pos = row['position']
            raw_score = row[score_col]
            if pd.notna(raw_score):
                try:
                    raw_score = float(raw_score)
                    if is_pvalue_col:
                        pval = max(raw_score, 1e-300)
                        score = -np.log10(pval)
                    else:
                        score = raw_score
                    pred_lookup[(tx, pos)] = score
                except (ValueError, TypeError):
                    pass

    # Build truth site lookup for efficient matching
    truth_pos_lookup = set()
    for _, row in truth_pos_subset.iterrows():
        truth_pos_lookup.add((row['transcript'], row['position']))

    truth_neg_lookup = set()
    if not truth_neg_subset.empty:
        for _, row in truth_neg_subset.iterrows():
            truth_neg_lookup.add((row['transcript'], row['position']))

    labels = []
    scores = []

    # For each covered site, determine label and score
    for tx, pos in all_covered_sites:
        # Check if this is a positive site (within window)
        is_positive = False
        for tpos in truth_pos_lookup:
            if tpos[0] == tx and abs(tpos[1] - pos) <= window:
                is_positive = True
                break

        # Check if this is a negative site (within window)
        is_negative = False
        if not is_positive:
            for tpos in truth_neg_lookup:
                if tpos[0] == tx and abs(tpos[1] - pos) <= window:
                    is_negative = True
                    break

        # Get score for this site
        if (tx, pos) in pred_lookup:
            score = pred_lookup[(tx, pos)]
        else:
            # Fair comparison: tool didn't call this site -> score=0 (negative prediction)
            score = 0.0

        if is_positive:
            labels.append(1)
            scores.append(score)
        elif is_negative or len(truth_neg_lookup) == 0:
            # Include as negative if explicit negative or if no explicit negatives
            # (in latter case, all non-positive predictions are treated as negatives)
            labels.append(0)
            scores.append(score)

    if len(set(labels)) < 2 or len(scores) == 0:
        return np.nan, np.nan

    try:
        from sklearn.metrics import roc_auc_score, average_precision_score
        auroc = roc_auc_score(labels, scores)
        auprc = average_precision_score(labels, scores)
        return auprc, auroc
    except Exception:
        return np.nan, np.nan


records = []
records_by_comparison = []  # For per-comparison metrics
records_by_negtype = {}  # tool -> list of per-negative-type metrics

# Get all modification types from positive sites
all_mod_types = set(truth_pos['modification_type'].unique())

# ── Helper: detect k-mer column ───────────────────────────────────────────────
def detect_kmer_column(df):
    """Detect the k-mer column in a dataframe."""
    kmer_cols = ['kmer', 'k_mer', 'k-mer', 'motif', 'sequence', 'context']
    for col in kmer_cols:
        if col in df.columns:
            return col
    # Case-insensitive search
    for col in df.columns:
        if col.lower() in ['kmer', 'k_mer', 'k-mer']:
            return col
    return None

# ── Helper: categorize negative sites ─────────────────────────────────────────
def categorize_negative_sites(pred_df, truth_pos_subset, kmer_col=None, window=0):
    """
    Categorize negative sites into 3 types based on their relationship to positive sites.

    Categories:
    1. All Other Sites: Positions not in positive sites (general negatives)
    2. Same K-mer Different Position: Positions sharing k-mer with positive sites
    3. Same Base Position: Positions at exact same location as positive sites
       (these are false positives at the modification site)

    Returns:
        dict: {neg_type: DataFrame of negative sites}
    """
    # Build positive site lookup
    truth_pos_set = set()
    truth_pos_kmers = {}  # (tx, pos) -> kmer
    for _, row in truth_pos_subset.iterrows():
        tx = row['transcript']
        pos = row['position']
        truth_pos_set.add((tx, pos))
        # Store k-mer if available
        if kmer_col and kmer_col in row:
            truth_pos_kmers[(tx, pos)] = row[kmer_col]

    # Get all kmers from positive sites for "same k-mer" category
    positive_kmers = set(truth_pos_kmers.values()) if truth_pos_kmers else set()

    neg_by_type = {
        'all_other': [],      # Type 1: All other sites (not positive)
        'same_kmer': [],      # Type 2: Same k-mer, different position
        'same_position': [],  # Type 3: Same base position (FP at exact pos)
    }

    for _, row in pred_df.iterrows():
        tx = row['transcript']
        pos = row['position']

        # Check if this is a positive site (within window)
        is_positive = False
        matched_pos = None
        for ttx, tpos in truth_pos_set:
            if ttx == tx and abs(tpos - pos) <= window:
                is_positive = True
                matched_pos = (ttx, tpos)
                break

        if is_positive:
            continue  # Skip positive sites

        # Categorize this negative
        pred_kmer = row.get(kmer_col) if kmer_col else None

        # Check if same position as a positive site (exact match but not positive)
        # This happens when tool calls at the exact position but it's not in truth
        is_same_position = False
        for ttx, tpos in truth_pos_set:
            if ttx == tx and tpos == pos:
                is_same_position = True
                break

        if is_same_position:
            neg_by_type['same_position'].append({
                'transcript': tx,
                'position': pos,
                'kmer': pred_kmer,
                'neg_type': 'same_position'
            })
        elif pred_kmer and pred_kmer in positive_kmers:
            neg_by_type['same_kmer'].append({
                'transcript': tx,
                'position': pos,
                'kmer': pred_kmer,
                'neg_type': 'same_kmer'
            })
        else:
            neg_by_type['all_other'].append({
                'transcript': tx,
                'position': pos,
                'kmer': pred_kmer,
                'neg_type': 'all_other'
            })

    # Convert to DataFrames
    result = {}
    for neg_type, sites in neg_by_type.items():
        if sites:
            result[neg_type] = pd.DataFrame(sites)
        else:
            result[neg_type] = pd.DataFrame(columns=['transcript', 'position', 'kmer', 'neg_type'])

    return result

# ── Build union of all covered sites per modification type (for fair comparison) ──
print("\nBuilding union of all covered sites for fair comparison...")
all_covered_sites_by_mod = {}  # mod_type -> set of (tx, pos) tuples

for mod_type in all_mod_types:
    all_covered_sites_by_mod[mod_type] = set()
    for tool, pred_df in tool_dfs.items():
        # Get relevant transcripts for this modification type
        truth_subset = truth_pos[truth_pos['modification_type'] == mod_type]
        relevant_tx = set(truth_subset['transcript'].unique())

        # Add sites from this tool
        tool_sites = pred_df[pred_df['transcript'].isin(relevant_tx)]
        for _, row in tool_sites.iterrows():
            all_covered_sites_by_mod[mod_type].add((row['transcript'], row['position']))

    print(f"  {mod_type}: {len(all_covered_sites_by_mod[mod_type])} unique sites covered by all tools")

for tool, pred_df in tool_dfs.items():
    for mod_type in all_mod_types:
        truth_subset = truth_pos[truth_pos['modification_type'] == mod_type].copy()

        if truth_subset.empty:
            continue

        # Get unique transcripts for this mod type
        all_tx = set(truth_subset['transcript'].unique())

        # Filter predictions to relevant transcripts
        pred_subset = pred_df[pred_df['transcript'].isin(all_tx)].copy()

        # Get all covered sites for this modification type
        all_covered_sites = all_covered_sites_by_mod.get(mod_type, set())

        # Auto-infer negative sites if not provided
        # For fair comparison, use the union of all covered sites as the universe
        truth_neg_subset = truth_neg[truth_neg['modification_type'] == mod_type].copy()
        is_inferred_negative = False

        if truth_neg_subset.empty and all_covered_sites:
            # Create inferred negative sites from all covered sites that don't match truth
            neg_sites = []
            truth_pos_set = set()
            for _, truth_row in truth_subset.iterrows():
                truth_pos_set.add((truth_row['transcript'], truth_row['position']))

            for tx, pos in all_covered_sites:
                # Check if this site matches any truth site
                is_match = False
                for ttx, tpos in truth_pos_set:
                    if ttx == tx and abs(tpos - pos) <= 0:
                        is_match = True
                        break
                if not is_match:
                    neg_sites.append({
                        'transcript': tx,
                        'position': pos,
                        'modification_type': mod_type,
                        'label': '-'
                    })
            if neg_sites:
                truth_neg_subset = pd.DataFrame(neg_sites)
                is_inferred_negative = True

        score_col = tool_score_cols.get(tool)

        for window in windows:
            # Track which predicted sites matched truth sites
            matched_pred_indices = set()

            # For each positive truth site, check if there's a prediction within window
            tp = 0
            for _, truth_row in truth_subset.iterrows():
                tx = truth_row['transcript']
                pos = truth_row['position']

                matches = pred_subset[
                    (pred_subset['transcript'] == tx) &
                    (pred_subset['position'].between(pos - window, pos + window))
                ]

                if not matches.empty:
                    tp += 1
                    matched_pred_indices.update(matches.index.tolist())

            fn = len(truth_subset) - tp

            # FAIR COMPARISON: FP includes sites this tool called that aren't truth sites
            # AND sites other tools called but this tool didn't (already counted as not in pred_subset)
            # The key change: total predictions = all covered sites, not just this tool's calls
            fp = len(pred_subset) - len(matched_pred_indices)

            called_sites = len(pred_subset)
            total_predicted = called_sites
            total_truth = len(truth_subset)

            # Calculate TN (true negatives)
            tn = 0
            if not truth_neg_subset.empty and not is_inferred_negative:
                for _, neg_row in truth_neg_subset.iterrows():
                    tx = neg_row['transcript']
                    pos = neg_row['position']

                    nearby_preds = pred_subset[
                        (pred_subset['transcript'] == tx) &
                        (pred_subset['position'].between(pos - window, pos + window))
                    ]

                    if nearby_preds.empty:
                        tn += 1

            total_negative = len(truth_neg_subset)

            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = (2 * precision * recall / (precision + recall)
                  if (precision + recall) > 0 else 0)

            if is_inferred_negative or total_negative == 0:
                specificity = np.nan
            else:
                specificity = tn / total_negative

            if is_inferred_negative:
                mcc = np.nan
            else:
                mcc = compute_mcc(tp, fp, fn, tn)

            # AUPRC and AUROC with fair comparison
            auprc, auroc = compute_ranking_metrics_with_fair_comparison(
                pred_subset, all_covered_sites, truth_subset, truth_neg_subset,
                score_col, window, tool_name=tool
            )

            records.append({
                "tool": tool,
                "modification_type": mod_type,
                "window": window,
                "precision": precision,
                "recall": recall,
                "f1": f1,
                "tp": tp, "fp": fp, "fn": fn, "tn": tn,
                "specificity": specificity,
                "mcc": mcc,
                "auprc": auprc,
                "auroc": auroc,
                "called_sites": called_sites,
                "total_truth": total_truth,
                "total_predicted": total_predicted,
                "total_negative": total_negative,
            })

        # ── Per-negative-type metrics (detailed benchmark by negative category) ──
        # When groundtruth only has positive sites, we categorize false positives into:
        #
        # NEGATIVE SITE CATEGORIES (for tools predicting modification sites):
        # ─────────────────────────────────────────────────────────────────────────
        # 1. all_other (Type 1): All other positions not in positive sites
        #    - These are general negative sites, the easiest category
        #    - Represents positions that have no relationship to true modification sites
        #
        # 2. same_kmer (Type 2): Positions sharing k-mer with positive sites but different position
        #    - These are "harder" negatives because they share the same sequence context
        #    - The k-mer is identical but the position differs
        #    - Tests tool's ability to distinguish exact modification position within the k-mer
        #
        # 3. same_position (Type 3): Positions at the exact same base position
        #    - These are false positives at the exact modification position
        #    - Tool called a site at the correct position but it wasn't in the truth set
        #    - Most confusing category - could indicate:
        #      a) False positives (tool is wrong)
        #      b) Missing truth sites (truth set is incomplete)
        # ─────────────────────────────────────────────────────────────────────────
        if tool not in records_by_negtype:
            records_by_negtype[tool] = []

        # Detect k-mer column for this tool
        kmer_col = detect_kmer_column(pred_subset)

        # Only compute per-negative-type if truth_neg is empty (no explicit negatives)
        # and we're using inferred negatives
        if truth_neg.empty and not truth_subset.empty:
            # Build positive site lookup (exact matches)
            truth_pos_set = set()
            for _, row in truth_subset.iterrows():
                truth_pos_set.add((row['transcript'], row['position']))

            # Build k-mer lookup from positive sites
            positive_kmers = set()
            positive_kmer_by_tx = {}  # (tx) -> set of kmers at positive sites
            if kmer_col and kmer_col in truth_subset.columns:
                for _, row in truth_subset.iterrows():
                    tx = row['transcript']
                    kmer = row.get(kmer_col)
                    if pd.notna(kmer):
                        positive_kmers.add(kmer)
                        if tx not in positive_kmer_by_tx:
                            positive_kmer_by_tx[tx] = set()
                        positive_kmer_by_tx[tx].add(kmer)

            # Categorize each predicted site that is NOT a true positive
            # These become our negative site categories
            fp_by_type = {
                'all_other': [],      # Type 1: General negatives
                'same_kmer': [],      # Type 2: Same k-mer, different position
                'same_position': [],  # Type 3: Same exact position (FP at truth site)
            }

            for idx, row in pred_subset.iterrows():
                tx = row['transcript']
                pos = row['position']

                # Check if this prediction is at the exact position of a positive site
                is_exact_match = (tx, pos) in truth_pos_set

                # If exact match, it's a true positive, skip for negative categorization
                if is_exact_match:
                    continue

                # This is a negative (either FP or TN depending on window)
                # Categorize by relationship to positive sites

                pred_kmer = row.get(kmer_col) if kmer_col else None

                # Type 3: Same position - this prediction is at a different transcript
                # or different position from positive sites, so not same_position
                # (We already excluded exact matches above)

                # Type 2: Same k-mer - check if this position shares k-mer with any positive site
                is_same_kmer = False
                if pred_kmer and pred_kmer in positive_kmers:
                    is_same_kmer = True

                if is_same_kmer:
                    fp_by_type['same_kmer'].append(idx)
                else:
                    fp_by_type['all_other'].append(idx)

            # Count negatives by type
            for neg_type in ['all_other', 'same_kmer', 'same_position']:
                neg_count = len(fp_by_type[neg_type])

                # Get subset of predictions for this negative type
                if fp_by_type[neg_type]:
                    pred_subset_negtype = pred_subset.loc[fp_by_type[neg_type]]
                else:
                    pred_subset_negtype = pd.DataFrame(columns=pred_subset.columns)

                for window in windows:
                    # For per-negative-type analysis:
                    # - TP = 0 (these are all negatives, no true positives in this subset)
                    # - FP = predictions in this subset that don't match any truth site
                    # - TN = all other negatives in this category not predicted

                    # Actually, we need to reconsider the approach:
                    # The user wants to see how the tool performs on different NEGATIVE SETS
                    # not different PREDICTION subsets.
                    #
                    # So we should compute:
                    # - For each negative type, treat those sites as the negative set
                    # - Compute precision/recall against that negative set
                    # - This shows tool's discrimination ability for each negative type

                    # Count FP in this negative type (predictions that are in this category)
                    fp_neg = len(pred_subset_negtype)

                    # Count TN = total negatives of this type - FP (negatives correctly not called)
                    tn_neg = neg_count - fp_neg if neg_count > 0 else 0

                    # For precision/recall, we use all truth positives and all tool predictions
                    # But specificity is computed against each negative type separately
                    precision_neg = precision  # Same as overall (uses all predictions)
                    recall_neg = recall  # Same as overall (uses all truth positives)
                    f1_neg = f1  # Same as overall

                    specificity_neg = tn_neg / neg_count if neg_count > 0 else np.nan

                    # AUPRC/AUROC for this negative type
                    # Create the negative set for this type
                    truth_neg_negtype = pd.DataFrame([
                        {'transcript': pred_subset.loc[idx, 'transcript'],
                         'position': pred_subset.loc[idx, 'position'],
                         'modification_type': mod_type, 'label': '-'}
                        for idx in fp_by_type[neg_type]
                    ]) if fp_by_type[neg_type] else pd.DataFrame()

                    # Build covered sites for fair comparison
                    covered_sites_neg = set()
                    for idx in fp_by_type[neg_type]:
                        r = pred_subset.loc[idx]
                        covered_sites_neg.add((r['transcript'], r['position']))

                    # Compute AUPRC/AUROC treating this negative type as the only negatives
                    auprc_neg, auroc_neg = compute_ranking_metrics_with_fair_comparison(
                        pred_subset, covered_sites_neg, truth_subset, truth_neg_negtype,
                        score_col, window, tool_name=tool
                    )

                    records_by_negtype[tool].append({
                        "tool": tool,
                        "modification_type": mod_type,
                        "negative_type": neg_type,
                        "window": window,
                        "precision": precision_neg,
                        "recall": recall_neg,
                        "f1": f1_neg,
                        "tp": tp,  # Same as overall
                        "fp": fp_neg,  # FP in this negative category
                        "fn": fn,  # Same as overall
                        "tn": tn_neg,  # TN in this negative category
                        "specificity": specificity_neg,
                        "auprc": auprc_neg,
                        "auroc": auroc_neg,
                        "called_sites": called_sites,
                        "total_truth": total_truth,
                        "total_negative": neg_count,
                        "total_predicted": total_predicted,
                        "fp_in_category": fp_neg,  # How many FPs fall into this category
                    })

        # ── Per-comparison metrics (for scatterplot) ──
        if tool in tool_comparisons:
            for comparison, comp_df in tool_comparisons[tool]:
                # Filter to relevant transcripts for this modification type
                comp_subset = comp_df[comp_df['transcript'].isin(all_tx)].copy()

                if comp_subset.empty:
                    continue

                # Get covered sites for this comparison
                comp_covered_sites = set()
                for _, row in comp_subset.iterrows():
                    comp_covered_sites.add((row['transcript'], row['position']))

                score_col = tool_score_cols.get(tool)

                for window in windows:
                    # Compute metrics for this comparison
                    matched_indices = set()
                    tp_comp = 0

                    for _, truth_row in truth_subset.iterrows():
                        tx = truth_row['transcript']
                        pos = truth_row['position']

                        matches = comp_subset[
                            (comp_subset['transcript'] == tx) &
                            (comp_subset['position'].between(pos - window, pos + window))
                        ]

                        if not matches.empty:
                            tp_comp += 1
                            matched_indices.update(matches.index.tolist())

                    fn_comp = len(truth_subset) - tp_comp
                    fp_comp = len(comp_subset) - len(matched_indices)

                    precision_comp = tp_comp / (tp_comp + fp_comp) if (tp_comp + fp_comp) > 0 else 0
                    recall_comp = tp_comp / (tp_comp + fn_comp) if (tp_comp + fn_comp) > 0 else 0
                    f1_comp = (2 * precision_comp * recall_comp / (precision_comp + recall_comp)
                               if (precision_comp + recall_comp) > 0 else 0)

                    # Calculate TN for this comparison
                    tn_comp = 0
                    if not truth_neg_subset.empty and not is_inferred_negative:
                        for _, neg_row in truth_neg_subset.iterrows():
                            tx = neg_row['transcript']
                            pos = neg_row['position']
                            nearby_preds = comp_subset[
                                (comp_subset['transcript'] == tx) &
                                (comp_subset['position'].between(pos - window, pos + window))
                            ]
                            if nearby_preds.empty:
                                tn_comp += 1

                    # Specificity and MCC
                    total_negative_comp = len(truth_neg_subset)
                    if is_inferred_negative or total_negative_comp == 0:
                        specificity_comp = np.nan
                        mcc_comp = np.nan
                    else:
                        specificity_comp = tn_comp / total_negative_comp
                        mcc_comp = compute_mcc(tp_comp, fp_comp, fn_comp, tn_comp)

                    # Compute AUPRC/AUROC for this comparison
                    auprc_comp, auroc_comp = compute_ranking_metrics_with_fair_comparison(
                        comp_subset, comp_covered_sites, truth_subset, truth_neg_subset,
                        score_col, window, tool_name=tool
                    )

                    records_by_comparison.append({
                        "tool": tool,
                        "modification_type": mod_type,
                        "comparison": comparison,
                        "window": window,
                        "precision": precision_comp,
                        "recall": recall_comp,
                        "f1": f1_comp,
                        "tp": tp_comp,
                        "fp": fp_comp,
                        "fn": fn_comp,
                        "tn": tn_comp,
                        "specificity": specificity_comp,
                        "mcc": mcc_comp,
                        "auprc": auprc_comp,
                        "auroc": auroc_comp,
                        "called_sites": len(comp_subset),
                        "total_truth": len(truth_subset),
                        "total_negative": total_negative_comp,
                    })

# Column definitions for output files
mod_type_columns = [
    "tool", "modification_type", "window", "precision", "recall",
    "f1", "tp", "fp", "fn", "tn", "specificity", "mcc",
    "auprc", "auroc", "called_sites",
    "total_truth", "total_predicted", "total_negative"
]
overall_columns = [
    "tool", "window", "precision", "recall",
    "f1", "tp", "fp", "fn", "tn", "specificity", "mcc",
    "auprc", "auroc", "called_sites",
    "total_truth", "total_predicted", "total_negative"
]
comparison_columns = [
    "tool", "modification_type", "comparison", "window", "precision", "recall",
    "f1", "tp", "fp", "fn", "tn", "specificity", "mcc", "auprc", "auroc",
    "called_sites", "total_truth", "total_negative"
]
negtype_columns = [
    "tool", "modification_type", "negative_type", "window", "precision", "recall",
    "f1", "tp", "fp", "fn", "tn", "specificity", "auprc", "auroc",
    "called_sites", "total_truth", "total_negative", "total_predicted"
]

# ── Write aggregated files (combining all tools) ───────────────────────────────────
print("\nWriting aggregated benchmark files...")

# Aggregate all per-tool records for each output type
all_records = []
all_comparison_records = []
all_negtype_records = []
all_overall_records = []

for tool in sorted(tool_dfs.keys()):
    # Collect per-modification-type records
    tool_records = [r for r in records if r["tool"] == tool]
    all_records.extend(tool_records)

    # Collect per-comparison records
    tool_comparison_records = [r for r in records_by_comparison if r["tool"] == tool]
    all_comparison_records.extend(tool_comparison_records)

    # Collect per-negative-type records
    tool_negtype_records = records_by_negtype.get(tool, [])
    all_negtype_records.extend(tool_negtype_records)

    # Compute overall records (aggregated across modification types for each tool)
    for window in windows:
        tool_window_records = [r for r in tool_records if r["window"] == window]
        if not tool_window_records:
            continue

        # Sum confusion matrix values
        tp_sum = sum(r["tp"] for r in tool_window_records)
        fp_sum = sum(r["fp"] for r in tool_window_records)
        fn_sum = sum(r["fn"] for r in tool_window_records)
        tn_sum = sum(r["tn"] for r in tool_window_records if not pd.isna(r["tn"]))
        tn_sum = tn_sum if tn_sum > 0 else 0

        total_truth_sum = sum(r["total_truth"] for r in tool_window_records)
        total_predicted_sum = sum(r["total_predicted"] for r in tool_window_records)
        total_negative_sum = sum(r["total_negative"] for r in tool_window_records)
        called_sites_sum = sum(r["called_sites"] for r in tool_window_records)

        # Check if any record has inferred negatives
        any_inferred = any(pd.isna(r.get("specificity")) for r in tool_window_records)

        # Calculate metrics
        precision = tp_sum / (tp_sum + fp_sum) if (tp_sum + fp_sum) > 0 else 0
        recall = tp_sum / (tp_sum + fn_sum) if (tp_sum + fn_sum) > 0 else 0
        f1 = (2 * precision * recall / (precision + recall)
              if (precision + recall) > 0 else 0)

        # Specificity and MCC
        if any_inferred or total_negative_sum == 0:
            specificity = np.nan
            mcc = np.nan
        else:
            specificity = tn_sum / total_negative_sum if total_negative_sum > 0 else np.nan
            mcc = compute_mcc(tp_sum, fp_sum, fn_sum, tn_sum)

        # AUPRC/AUROC: weighted average across modification types
        auprc_vals = [r["auprc"] for r in tool_window_records if not pd.isna(r["auprc"])]
        auroc_vals = [r["auroc"] for r in tool_window_records if not pd.isna(r["auroc"])]
        auprc_avg = np.mean(auprc_vals) if auprc_vals else np.nan
        auroc_avg = np.mean(auroc_vals) if auroc_vals else np.nan

        all_overall_records.append({
            "tool": tool,
            "window": window,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "tp": tp_sum, "fp": fp_sum, "fn": fn_sum, "tn": tn_sum,
            "specificity": specificity,
            "mcc": mcc,
            "auprc": auprc_avg,
            "auroc": auroc_avg,
            "called_sites": called_sites_sum,
            "total_truth": total_truth_sum,
            "total_predicted": total_predicted_sum,
            "total_negative": total_negative_sum,
        })

# Write aggregated accuracy_summary.tsv (by modification type)
if "by_mod_type" in aggregated_files:
    if all_records:
        agg_df = pd.DataFrame(all_records).sort_values(
            ["tool", "modification_type", "window", "f1"],
            ascending=[True, True, True, False]
        )
    else:
        agg_df = pd.DataFrame(columns=mod_type_columns)
    os.makedirs(os.path.dirname(aggregated_files["by_mod_type"]), exist_ok=True)
    agg_df.to_csv(aggregated_files["by_mod_type"], sep='\t', index=False)
    print(f"Saved aggregated per-modification-type metrics to {aggregated_files['by_mod_type']}")

# Write aggregated accuracy_summary_overall.tsv
if "overall" in aggregated_files:
    if all_overall_records:
        overall_df = pd.DataFrame(all_overall_records).sort_values(
            ["tool", "window", "f1"],
            ascending=[True, True, False]
        )
    else:
        overall_df = pd.DataFrame(columns=overall_columns)
    os.makedirs(os.path.dirname(aggregated_files["overall"]), exist_ok=True)
    overall_df.to_csv(aggregated_files["overall"], sep='\t', index=False)
    print(f"Saved aggregated overall metrics to {aggregated_files['overall']}")

# Write aggregated accuracy_summary_by_comparison.tsv
if "by_comparison" in aggregated_files:
    if all_comparison_records:
        comp_df = pd.DataFrame(all_comparison_records).sort_values(
            ["tool", "comparison", "window", "f1"],
            ascending=[True, True, True, False]
        )
    else:
        comp_df = pd.DataFrame(columns=comparison_columns)
    os.makedirs(os.path.dirname(aggregated_files["by_comparison"]), exist_ok=True)
    comp_df.to_csv(aggregated_files["by_comparison"], sep='\t', index=False)
    print(f"Saved aggregated per-comparison metrics to {aggregated_files['by_comparison']}")

# Write aggregated accuracy_summary_by_negative_type.tsv
if "by_negative_type" in aggregated_files:
    if all_negtype_records:
        negtype_df = pd.DataFrame(all_negtype_records).sort_values(
            ["tool", "modification_type", "negative_type", "window", "f1"],
            ascending=[True, True, True, True, False]
        )
    else:
        negtype_df = pd.DataFrame(columns=negtype_columns)
    os.makedirs(os.path.dirname(aggregated_files["by_negative_type"]), exist_ok=True)
    negtype_df.to_csv(aggregated_files["by_negative_type"], sep='\t', index=False)
    print(f"Saved aggregated per-negative-type metrics to {aggregated_files['by_negative_type']}")

# ── Write dedicated count tables ─────────────────────────────────────────────
# Create focused count tables for easier analysis of site calling behavior
benchmark_dir = os.path.dirname(aggregated_files.get("by_mod_type", ""))

# 1. called_sites_by_comparison.tsv - site counts per tool per comparison
if all_comparison_records:
    count_cols = ["tool", "modification_type", "comparison", "window",
                  "tp", "fp", "fn", "tn", "called_sites", "total_truth", "total_negative"]
    count_df = pd.DataFrame(all_comparison_records)[count_cols].sort_values(
        ["tool", "modification_type", "comparison", "window"]
    )
    count_path = os.path.join(benchmark_dir, "called_sites_by_comparison.tsv")
    count_df.to_csv(count_path, sep='\t', index=False)
    print(f"Saved called sites by comparison to {count_path}")

# 2. called_sites_summary.tsv - aggregated site counts per tool
if all_overall_records:
    summary_count_cols = ["tool", "window", "tp", "fp", "fn", "tn",
                          "called_sites", "total_truth", "total_predicted", "total_negative"]
    summary_count_df = pd.DataFrame(all_overall_records)[summary_count_cols].sort_values(
        ["tool", "window"]
    )
    summary_count_path = os.path.join(benchmark_dir, "called_sites_summary.tsv")
    summary_count_df.to_csv(summary_count_path, sep='\t', index=False)
    print(f"Saved called sites summary to {summary_count_path}")

print("\nBenchmark accuracy analysis complete!")
