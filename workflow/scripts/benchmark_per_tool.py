import os
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

# Import shared utilities
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_utils import tool_from_path, normalize_columns, is_pvalue_column, match_column


TOOL_SCORE_MAP = {
    'xpore': ['pval_CASE_vs_CONTROL', 'z_score_CASE_vs_CONTROL', 'diff_mod_rate_CASE_vs_CONTROL',
              'p_value', 'adjusted_p_value', 'statistic', 'pval', 'padj'],
    'nanocompore': ['GMM_logit_pvalue', 'KS_dwell_pvalue', 'KS_intensity_pvalue', 'Logit_LOR',
                    'p_value_glm', 'p_value_ks', 'p_value'],
    'baleen': ['p_value', 'padj', 'z_test', 'mod_ratio_diff', 'mod_ratio_fc', 'nLOD_diff',
               'native_mod_ratio', 'native_normalized_LOD', 'native_LOD',
               'pvalue', 'effect_size', 'score'],
    'pybaleen': ['pvalue', 'padj', 'effect_size', 'p_value', 'adjusted_p_value', 'stoichiometry'],
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


def evaluate_column(scores_arr, labels, col_name, is_pval):
    if is_pval:
        scores_arr = -np.log10(np.clip(scores_arr, 1e-300, None))
        transform = '-log10'
    else:
        transform = 'none'

    if len(np.unique(labels)) < 2:
        return {
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
        }

    try:
        auroc = roc_auc_score(labels, scores_arr)
    except Exception:
        auroc = np.nan

    try:
        prauc = average_precision_score(labels, scores_arr)
    except Exception:
        prauc = np.nan

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

    if is_pval:
        best_thresh_original = float(10 ** (-best_thresh))
    else:
        best_thresh_original = float(best_thresh)

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
        # Fair mode: uncalled sites get the minimum possible score so they
        # are always classified as negative regardless of threshold.
        # This correctly counts them as TN (truth=0) or FN (truth=1).
        valid_scores = scores_arr[valid_mask]
        if len(valid_scores) > 0:
            min_score = np.min(valid_scores) - 1.0
        else:
            min_score = -1.0
        scores_arr = np.where(np.isnan(scores_arr), min_score, scores_arr)
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
