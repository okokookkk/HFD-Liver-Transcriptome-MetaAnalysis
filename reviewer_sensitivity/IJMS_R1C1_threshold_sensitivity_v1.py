#!/usr/bin/env python3
"""IJMS Reviewer 1, Major Comment 1: three-criterion threshold sensitivity.

Runs on the existing mouse gene-level multilevel meta-analysis output.
No DESeq2 or metafor model is refitted.

Example:
    python3 IJMS_R1C1_threshold_sensitivity_v1.py \
        --root .

Outputs (default: 00_metadata/analysis_outputs/revision_R1C1_thresholds):
    threshold_one_at_a_time_v1.tsv
    threshold_grid_27_v1.tsv
    original_core_59_stringency_v1.tsv
    audit_summary_v1.txt
"""
import argparse
import csv
import gzip
import math
from itertools import product
from pathlib import Path

# Z quantiles used by the original R analysis: qnorm((1 + level)/2)
Z_CRITICAL = {0.90: 1.6448536269514722,
              0.95: 1.959963984540054,
              0.99: 2.5758293035489004}
FDR_LEVELS = (0.01, 0.05, 0.10)
PI_LEVELS = (0.90, 0.95, 0.99)
CONS_LEVELS = (0.75, 0.80, 0.85)
CORE_BASE = (0.05, 0.95, 0.80)


def open_text(path):
    return gzip.open(path, 'rt', newline='') if path.name.endswith('.gz') else path.open('r', newline='')


def read_table(path):
    with open_text(path) as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def save_tsv(path, headers, rows):
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def require_float(record, key, gene):
    try:
        value = float(record[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f'{gene}: invalid {key}: {record.get(key)!r}') from exc
    if not math.isfinite(value):
        raise ValueError(f'{gene}: non-finite {key}: {value!r}')
    return value


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, default=Path('.'), help='HFD_MetaAnalysis directory')
    ap.add_argument('--outdir', type=Path, default=None, help='output directory')
    args = ap.parse_args()
    root = args.root.resolve()
    outdir = (args.outdir or root / '00_metadata/analysis_outputs/revision_R1C1_thresholds').resolve()

    meta_path = root / '05_meta/gene_level_v1_2/HFD_gene_multilevel_meta_primary_v1_2.tsv.gz'
    core_path = root / '05_meta/stable_core_v1/HFD_stable_core_final_v1.tsv'
    for path in (meta_path, core_path):
        if not path.is_file():
            ap.error(f'File not found: {path}')

    raw = read_table(meta_path)
    original_core = read_table(core_path)
    core_ids = {r['Geneid'] for r in original_core}
    if len(core_ids) != len(original_core):
        raise ValueError('Duplicate Geneid in original core table')
    if len({r['Geneid'] for r in raw}) != len(raw):
        raise ValueError('Duplicate Geneid in meta-analysis table')

    records = []
    stored_pi_disagreement = []
    bad_components = []
    for r in raw:
        g = r['Geneid']
        effect = require_float(r, 'pooled_log2FC', g)
        se = require_float(r, 'pooled_SE', g)
        fdr = require_float(r, 'pooled_FDR_BH', g)
        tau_gse = require_float(r, 'tau2_GSE', g)
        tau_contrast = require_float(r, 'tau2_contrast_within_GSE', g)
        cons = require_float(r, 'GSE_balanced_direction_consistency', g)
        if se < 0 or tau_gse < -1e-12 or tau_contrast < -1e-12 or not (0 <= cons <= 1) or not (0 <= fdr <= 1):
            bad_components.append(g)
            continue
        total_variance = se * se + max(0.0, tau_gse) + max(0.0, tau_contrast)
        pred_sd = math.sqrt(total_variance)
        pi_pass = {p: ((effect - z * pred_sd) > 0 or (effect + z * pred_sd) < 0)
                   for p, z in Z_CRITICAL.items()}
        recorded_pass = str(r.get('prediction_interval_excludes_zero', '')).strip().upper() == 'YES'
        if pi_pass[0.95] != recorded_pass:
            stored_pi_disagreement.append(g)
        records.append({'Geneid': g, 'fdr': fdr, 'cons': cons, 'pred_sd': pred_sd,
                        'effect': effect, 'pi_pass': pi_pass,
                        'N_GSE': r.get('N_independent_GSE', '')})
    if bad_components:
        raise ValueError(f'Unexpected negative variance or invalid scaling in {len(bad_components)} genes. First: {bad_components[:5]}')
    if stored_pi_disagreement:
        raise ValueError('Recomputed 95% PI exclusion disagrees with original for '
                         f'{len(stored_pi_disagreement)} genes: {stored_pi_disagreement[:10]}')

    def select(fdr_thr, pi_thr, cons_thr):
        return {x['Geneid'] for x in records
                if x['fdr'] < fdr_thr and x['pi_pass'][pi_thr] and x['cons'] >= cons_thr}

    original = select(*CORE_BASE)
    if len(original) != 59 or original != core_ids:
        missing = sorted(core_ids - original)
        extra = sorted(original - core_ids)
        raise ValueError(f'Primary three-criterion core not reproduced: derived={len(original)}, '
                         f'final={len(core_ids)}, final-only={missing[:10]}, derived-only={extra[:10]}')

    def summary_row(kind, fdr, pi, cons):
        selected = select(fdr, pi, cons)
        overlap = len(selected & original)
        union = len(selected | original)
        return {'Analysis': kind, 'FDR_cutoff': fdr, 'PI_level': pi,
                'Directional_consistency_cutoff': cons,
                'Selected_genes': len(selected),
                'Overlap_with_original_59': overlap,
                'Original_59_retained_percent': round(overlap / len(original) * 100, 2),
                'New_genes_vs_original_59': len(selected - original),
                'Jaccard_vs_original_59': round(overlap / union, 4) if union else 1.0}

    axis = []
    for v in FDR_LEVELS:
        axis.append(summary_row('Vary_FDR_only', v, 0.95, 0.80))
    for v in PI_LEVELS:
        axis.append(summary_row('Vary_PI_only', 0.05, v, 0.80))
    for v in (0.70, 0.75, 0.80, 0.85, 0.90):
        axis.append(summary_row('Vary_consistency_only', 0.05, 0.95, v))

    grid = [summary_row('Three_factor_grid', f, p, c)
            for f, p, c in product(FDR_LEVELS, PI_LEVELS, CONS_LEVELS)]
    columns = ['Analysis', 'FDR_cutoff', 'PI_level', 'Directional_consistency_cutoff',
               'Selected_genes', 'Overlap_with_original_59',
               'Original_59_retained_percent', 'New_genes_vs_original_59',
               'Jaccard_vs_original_59']
    outdir.mkdir(parents=True, exist_ok=True)
    save_tsv(outdir / 'threshold_one_at_a_time_v1.tsv', columns, axis)
    save_tsv(outdir / 'threshold_grid_27_v1.tsv', columns, grid)

    strict_sets = {
        'FDR_0.01': select(0.01, 0.95, 0.80),
        'PI_99_percent': select(0.05, 0.99, 0.80),
        'Consistency_0.85': select(0.05, 0.95, 0.85),
        'All_three_stricter': select(0.01, 0.99, 0.85),
    }
    by_gene = {x['Geneid']: x for x in records}
    core_rows = []
    for gene in sorted(original):
        x = by_gene[gene]
        core_rows.append({
            'Geneid': gene, 'N_independent_GSE': x['N_GSE'],
            'Pooled_FDR_BH': x['fdr'],
            'GSE_balanced_consistency': x['cons'],
            **{f'Retained_under_{key}': int(gene in value)
               for key, value in strict_sets.items()}
        })
    save_tsv(outdir / 'original_core_59_stringency_v1.tsv',
             ['Geneid', 'N_independent_GSE', 'Pooled_FDR_BH', 'GSE_balanced_consistency'] +
             [f'Retained_under_{key}' for key in strict_sets], core_rows)

    notes = [
        'IJMS Reviewer 1 Major Comment 1: threshold sensitivity audit',
        f'Meta-analysis genes={len(records)}',
        f'Original core genes={len(original)}',
        'Original three-criterion definition: FDR<0.05; normal-approximation 95% PI excludes 0; GSE-balanced consistency>=0.80.',
        'Original core IDs exactly match HFD_stable_core_final_v1.tsv: YES',
        'All stored 95% PI exclusion flags agree with recomputed PI using pooled_SE^2 + tau2_GSE + tau2_contrast_within_GSE: YES',
        'These are threshold reclassifications based on existing fixed fitted models; no new meta-analysis model was fitted.',
        'A changed number of selected genes does not alone demonstrate biological robustness.',
        'The 27-condition grid is descriptive/exploratory, not a search for an optimal new cutoff.',
    ]
    (outdir / 'audit_summary_v1.txt').write_text('\n'.join(notes) + '\n', encoding='utf-8')
    print('\n'.join(notes))
    print('\nONE-AT-A-TIME SENSITIVITY:')
    print('Axis\tFDR\tPI\tConsistency\tN_selected\tOf_original_59_retained\tNew_vs_59')
    for row in axis:
        print(f"{row['Analysis']}\t{row['FDR_cutoff']:.2f}\t{row['PI_level']:.2f}\t"
              f"{row['Directional_consistency_cutoff']:.2f}\t{row['Selected_genes']}\t"
              f"{row['Overlap_with_original_59']}\t{row['New_genes_vs_original_59']}")
    print('\nSTRICTER SINGLE CRITERIA AND COMBINED:')
    for key, sel in strict_sets.items():
        print(f'{key}: {len(sel)} of original 59 remain')
    print(f'\nOutputs: {outdir}')


if __name__ == '__main__':
    main()
