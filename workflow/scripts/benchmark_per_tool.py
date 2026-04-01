import os
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve



def tool_from_path(path):
    """Extract tool name from file path."""
    parts = path.split(os.sep)
    for i, part in enumerate(parts):
        if part == "modifications" and i + 1 < len(parts):
            return parts[i + 1]
    return os.path.basename(os.path.dirname(os.path.dirname(path)))


def normalize_columns(df):
    """Normalize various column names to standard 'transcript' and 'position'."""
    col_mapping = {}
    transcript_cols = ['transcript', 'id', 'ref_id', 'chrom', 'transcript_id']
    for col in transcript_cols:
        if col in df.columns:
            col_mapping[col] = 'transcript'
            break
    position_cols = ['position', 'pos', 'start', 'start_loc', 'transcript_pos', 'transcript_loc']
    for col in position_cols:
        if col in df.columns and col not in col_mapping.values():
            col_mapping[col] = 'position'
            break
    if col_mapping:
        df = df.rename(columns=col_mapping)
    return df


def match_column(df_columns, pattern):
    """Match a pattern against DataFrame columns using multiple strategies."""
    pattern_lower = pattern.lower().replace(' ', '_')
    if pattern in df_columns:
        return pattern
    for col in df_columns:
        if col.lower().replace(' ', '_') == pattern_lower:
            return col
    for col in df_columns:
        col_lower = col.lower().replace(' ', '_')
        if pattern_lower in col_lower:
            return col
    for col in df_columns:
        col_lower = col.lower().replace(' ', '_')
        if col_lower in pattern_lower:
            return col
    return None


def is_pvalue_column(col_name):
    """Determine if a column is a p-value (lower=better) or score (higher=better)."""
    if not col_name:
        return False
    col_lower = col_name.lower().replace(' ', '_')
    if '-log10' in col_lower:
        return False
    if any(x in col_lower for x in ['p_value', 'pvalue', 'pval', 'padj', 'fdr', 'q_value', 'qval']):
        return True
    if any(x in col_lower for x in ['probability', 'prob', 'stoichiometry', 'logit',
                                      'z_score', 'score', 'delta', 'diff_mod', 'mod_ratio']):
        return False
    return False


TOOL_SCORE_MAP = {
    'xpore': ['pval_CASE_vs_CONTROL', 'z_score_CASE_vs_CONTROL', 'diff_mod_rate_CASE_vs_CONTROL',
              'p_value', 'adjusted_p_value', 'statistic', 'pval', 'padj'],
    'nanocompore': ['GMM_logit_pvalue', 'KS_dwell_pvalue', 'KS_intensity_pvalue', 'Logit_LOR',
                    'p_value_glm', 'p_value_ks', 'p_value'],
    'baleen': ['p_value', 'padj', 'z_test', 'mod_ratio_diff', 'mod_ratio_fc', 'nLOD_diff',
               'native_mod_ratio', 'native_normalized_LOD', 'native_LOD',
               'pvalue', 'effect_size', 'score'],
    'pybaleen': ['mod_ratio', 'pvalue', 'padj', 'effect_size', 'mean_p_mod',
                  'p_value', 'adjusted_p_value', 'stoichiometry'],
    'differr': ['-log10 P value', '-log10 FDR', 'odds ratio', 'G statistic', 'score'],
    'drummer': ['OR_padj', 'G_padj', 'G_test', 'log2_(OR)', 'p_value', 'pvalue', 'padj', 'odds_ratio'],
    'eligos2': ['pval', 'adjPval', 'p_value', 'q_value', 'error_rate', 'pvalue', 'qvalue'],
    'epinano': ['z_scores', 'delta_sum_err', 'z_score_prediction', 'score'],
    'psipore': ['p_value', 'psi_score', 'pvalue', 'score'],
}

GENERIC_SCORE_PATTERNS = [
    'p_value', 'pvalue', 'pval', '-log10_pvalue', '-log10 P value',
    'score', 'probability', 'prob', 'mod_score', 'mod_prob',
    'logit_pvalue', 'logit', 'z_score', 'z_scores',
]

OUTPUT_COLUMNS = [
    'score_column', 'is_pvalue', 'transform', 'auroc', 'prauc',
    'best_threshold', 'best_threshold_original', 'f1', 'precision', 'recall',
    'youden_threshold', 'youden_threshold_original', 'youden_j',
    'youden_sensitivity', 'youden_specificity',
    'n_sites', 'n_positive', 'n_negative',
]


def get_candidate_columns(df, tool_name):
    patterns = TOOL_SCORE_MAP.get(tool_name, GENERIC_SCORE_PATTERNS)
    df_columns = list(df.columns)
    seen = []
    result = []
    for pattern in patterns:
        matched = match_column(df_columns, pattern)
        if matched and matched not in seen:
            seen.append(matched)
            result.append((matched, is_pvalue_column(matched)))
    return result


def build_labels(eval_sites, truth_pos, max_window):
    truth_dict = {}
    for _, row in truth_pos.iterrows():
        tx = row['transcript']
        pos = int(row['position'])
        if tx not in truth_dict:
            truth_dict[tx] = []
        truth_dict[tx].append(pos)

    labels = []
    for tx, pos in eval_sites:
        pos = int(pos)
        label = 0
        if tx in truth_dict:
            for tp in truth_dict[tx]:
                if abs(pos - tp) <= max_window:
                    label = 1
                    break
        labels.append(label)
    return np.array(labels)


def compute_f1_at_threshold(scores, labels, threshold):
    preds = (scores >= threshold).astype(int)
    tp = int(np.sum((preds == 1) & (labels == 1)))
    fp = int(np.sum((preds == 1) & (labels == 0)))
    fn = int(np.sum((preds == 0) & (labels == 1)))
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return f1, precision, recall


INF_CAP = 1000.0  # Replace inf/-inf with a large finite value


def evaluate_column(scores_arr, labels, col_name, is_pval):
    if is_pval:
        scores_arr = -np.log10(np.clip(scores_arr, 1e-300, None))
        transform = '-log10'
    else:
        transform = 'none'

    # Cap inf values (e.g. differr's "-log10 P value" column contains inf from -log10(0))
    scores_arr = np.where(np.isposinf(scores_arr), INF_CAP, scores_arr)
    scores_arr = np.where(np.isneginf(scores_arr), -INF_CAP, scores_arr)

    n_sites = len(labels)
    n_positive = int(np.sum(labels == 1))
    n_negative = int(np.sum(labels == 0))

    empty_result = {
        'score_column': col_name,
        'is_pvalue': is_pval,
        'transform': transform,
        'auroc': np.nan,
        'prauc': np.nan,
        'best_threshold': np.nan,
        'best_threshold_original': np.nan,
        'f1': 0.0,
        'precision': 0.0,
        'recall': 0.0,
        'youden_threshold': np.nan,
        'youden_threshold_original': np.nan,
        'youden_j': np.nan,
        'youden_sensitivity': np.nan,
        'youden_specificity': np.nan,
        'n_sites': n_sites,
        'n_positive': n_positive,
        'n_negative': n_negative,
    }

    if len(np.unique(labels)) < 2:
        return empty_result

    try:
        auroc = roc_auc_score(labels, scores_arr)
    except Exception:
        auroc = np.nan

    try:
        prauc = average_precision_score(labels, scores_arr)
    except Exception:
        prauc = np.nan

    # --- F1-optimal threshold (from grid search) ---
    thresholds = np.linspace(scores_arr.min(), scores_arr.max(), 100)
    best_f1 = -1.0
    best_thresh = thresholds[0]
    best_precision = 0.0
    best_recall = 0.0

    for thresh in thresholds:
        f1, precision, recall = compute_f1_at_threshold(scores_arr, labels, thresh)
        if f1 > best_f1:
            best_f1 = f1
            best_thresh = thresh
            best_precision = precision
            best_recall = recall

    # Map threshold back to original score space
    if is_pval:
        best_thresh_original = float(10 ** (-best_thresh))
    else:
        best_thresh_original = float(best_thresh)

    # --- Youden's J threshold (from ROC curve) ---
    # J = sensitivity + specificity - 1 = TPR - FPR
    youden_thresh = np.nan
    youden_thresh_original = np.nan
    youden_j = np.nan
    youden_sensitivity = np.nan
    youden_specificity = np.nan

    try:
        fpr, tpr, roc_thresholds = roc_curve(labels, scores_arr)
        j_scores = tpr - fpr
        best_j_idx = np.argmax(j_scores)
        youden_j = float(j_scores[best_j_idx])
        youden_sensitivity = float(tpr[best_j_idx])
        youden_specificity = float(1 - fpr[best_j_idx])
        youden_thresh = float(roc_thresholds[best_j_idx])
        if is_pval:
            youden_thresh_original = float(10 ** (-youden_thresh))
        else:
            youden_thresh_original = float(youden_thresh)
    except Exception:
        pass

    return {
        'score_column': col_name,
        'is_pvalue': is_pval,
        'transform': transform,
        'auroc': auroc,
        'prauc': prauc,
        'best_threshold': float(best_thresh),
        'best_threshold_original': best_thresh_original,
        'f1': best_f1,
        'precision': best_precision,
        'recall': best_recall,
        'youden_threshold': youden_thresh,
        'youden_threshold_original': youden_thresh_original,
        'youden_j': youden_j,
        'youden_sensitivity': youden_sensitivity,
        'youden_specificity': youden_specificity,
        'n_sites': n_sites,
        'n_positive': n_positive,
        'n_negative': n_negative,
    }


def write_empty(output_scores, output_best_metrics, output_best_score):
    empty_df = pd.DataFrame(columns=OUTPUT_COLUMNS)
    empty_df.to_csv(output_scores, sep='\t', index=False)
    empty_df.to_csv(output_best_metrics, sep='\t', index=False)
    with open(output_best_score, 'w') as fh:
        fh.write('')


window_param = snakemake.params.window
if isinstance(window_param, list):
    max_window = max(int(w) for w in window_param)
else:
    max_window = int(window_param)

results_path = snakemake.input.results
truth_path = snakemake.input.truth_set
output_scores = snakemake.output.scores
output_best_metrics = snakemake.output.best_metrics
output_best_score = snakemake.output.best_score

fair_mode = 'union' in snakemake.input.keys()
has_negatives = 'negatives' in snakemake.input.keys()

tool_name = tool_from_path(results_path)

try:
    tool_df = pd.read_csv(results_path, sep='\t')
    tool_df = normalize_columns(tool_df)
except Exception:
    write_empty(output_scores, output_best_metrics, output_best_score)
    raise SystemExit(0)

if 'transcript' not in tool_df.columns or 'position' not in tool_df.columns:
    write_empty(output_scores, output_best_metrics, output_best_score)
    raise SystemExit(0)

tool_df['position'] = pd.to_numeric(tool_df['position'], errors='coerce')
tool_df = tool_df.dropna(subset=['position'])
tool_df['position'] = tool_df['position'].astype(int)

truth_df = pd.read_csv(truth_path, sep='\t')
truth_df = normalize_columns(truth_df)

if 'label' in truth_df.columns:
    truth_pos = truth_df[truth_df['label'] != '-'].copy()
else:
    truth_pos = truth_df.copy()

if truth_pos.empty:
    write_empty(output_scores, output_best_metrics, output_best_score)
    raise SystemExit(0)

truth_pos['position'] = pd.to_numeric(truth_pos['position'], errors='coerce')
truth_pos = truth_pos.dropna(subset=['position'])
truth_pos['position'] = truth_pos['position'].astype(int)

if fair_mode:
    union_df = pd.read_csv(snakemake.input.union, sep='\t')
    union_df = normalize_columns(union_df)
    union_df['position'] = pd.to_numeric(union_df['position'], errors='coerce')
    union_df = union_df.dropna(subset=['position'])
    union_df['position'] = union_df['position'].astype(int)
    eval_sites = list(zip(union_df['transcript'], union_df['position']))
else:
    truth_transcripts = set(truth_pos['transcript'].unique())
    tool_on_truth = tool_df[tool_df['transcript'].isin(truth_transcripts)]
    eval_sites = list(zip(tool_on_truth['transcript'], tool_on_truth['position']))

if has_negatives:
    neg_df = pd.read_csv(snakemake.input.negatives, sep='\t')
    neg_df = normalize_columns(neg_df)
    neg_df['position'] = pd.to_numeric(neg_df['position'], errors='coerce')
    neg_df = neg_df.dropna(subset=['position'])
    neg_df['position'] = neg_df['position'].astype(int)
    neg_sites = set(zip(neg_df['transcript'], neg_df['position']))

    # Build per-transcript positive positions for fast window lookup
    pos_by_tx = {}
    for _, row in truth_pos.iterrows():
        pos_by_tx.setdefault(row['transcript'], []).append(int(row['position']))

    # Only keep eval sites that are in truth positives (within window) or negatives
    filtered = []
    for tx, pos in eval_sites:
        pos_int = int(pos)
        is_pos = any(abs(pos_int - tp) <= max_window for tp in pos_by_tx.get(tx, []))
        if is_pos or (tx, pos_int) in neg_sites:
            filtered.append((tx, pos))
    eval_sites = filtered

if not eval_sites:
    write_empty(output_scores, output_best_metrics, output_best_score)
    raise SystemExit(0)

labels = build_labels(eval_sites, truth_pos, max_window)

tool_index = {(row['transcript'], row['position']): i
              for i, (_, row) in enumerate(tool_df.iterrows())}

candidate_cols = get_candidate_columns(tool_df, tool_name)

if not candidate_cols:
    write_empty(output_scores, output_best_metrics, output_best_score)
    raise SystemExit(0)

all_results = []

for col_name, is_pval in candidate_cols:
    raw_scores = []
    has_prediction = []
    for tx, pos in eval_sites:
        key = (tx, int(pos))
        if key in tool_index:
            row_idx = tool_index[key]
            val = tool_df.iloc[row_idx][col_name]
            try:
                val = float(val)
            except (ValueError, TypeError):
                val = np.nan
            raw_scores.append(val)
            has_prediction.append(True)
        else:
            raw_scores.append(np.nan)
            has_prediction.append(False)

    scores_arr = np.array(raw_scores, dtype=float)
    has_pred_arr = np.array(has_prediction)
    valid_mask = ~np.isnan(scores_arr)

    if not fair_mode:
        scores_arr = scores_arr[valid_mask]
        labels_eval = labels[valid_mask]
    else:
        # Fair mode: uncalled sites must get the worst possible score so they
        # are always classified as negative regardless of threshold.
        # For p-values (lower = more significant): fill with 1.0 (not significant),
        #   then -log10(1.0) = 0 (lowest after transform, no NaN/Inf issues).
        # For regular scores (higher = more modified): fill with min - 1.0
        #   so they always rank below any called site.
        if is_pval:
            fill_value = 1.0
        else:
            valid_scores = scores_arr[valid_mask]
            if len(valid_scores) > 0:
                fill_value = np.min(valid_scores) - 1.0
            else:
                fill_value = -1.0
        scores_arr = np.where(np.isnan(scores_arr), fill_value, scores_arr)
        labels_eval = labels

    if len(scores_arr) == 0:
        continue

    result = evaluate_column(scores_arr, labels_eval, col_name, is_pval)
    all_results.append(result)

if not all_results:
    write_empty(output_scores, output_best_metrics, output_best_score)
    raise SystemExit(0)

results_df = pd.DataFrame(all_results, columns=OUTPUT_COLUMNS)
results_df = results_df.sort_values('auroc', ascending=False).reset_index(drop=True)

results_df.to_csv(output_scores, sep='\t', index=False)

best_row = results_df.iloc[[0]]
best_row.to_csv(output_best_metrics, sep='\t', index=False)

best_col = results_df.iloc[0]['score_column']
with open(output_best_score, 'w') as fh:
    fh.write(str(best_col))
