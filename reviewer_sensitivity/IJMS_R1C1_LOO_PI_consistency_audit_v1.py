#!/usr/bin/env python3
"""Audit leave-one-contributing-GSE-out stability of the original 59-gene core.

Reuses existing rma.mv LOO refits; does NOT refit models, recompute genome-wide
BH-FDR, or establish complete three-criterion LOO membership.

Usage:
  python3 IJMS_R1C1_LOO_PI_consistency_audit_v1.py --root .

Outputs:
  00_metadata/analysis_outputs/revision_R1C1_LOO_criteria/
      loo_pi_consistency_detail_v1.tsv
      loo_pi_consistency_by_gene_v1.tsv
      audit_summary_v1.txt
"""
import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

Z95 = 1.959963984540054  # R qnorm(0.975)


def read_tsv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', newline='') as f:
        yield from csv.DictReader(f, delimiter='\t')


def write_tsv(path, rows, columns):
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def number(row, key):
    try:
        x = float(row[key])
    except (KeyError, ValueError, TypeError) as e:
        raise ValueError(f"Invalid numeric field {key} in row: {row}") from e
    if not math.isfinite(x):
        raise ValueError(f"Non-finite value in {key} for {row.get('Geneid')}")
    return x


def sign(x):
    return (x > 0) - (x < 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('.'))
    parser.add_argument('--outdir', type=Path, default=None)
    args = parser.parse_args()
    root = args.root.resolve()
    core_path = root / '05_meta/stable_core_v1/HFD_stable_core_final_v1.tsv'
    loo_path = root / '05_meta/stable_core_v1/HFD_stable_core_preLOO_leave_one_GSE_out_detail_v1.tsv.gz'
    effects_path = root / '00_metadata/analysis_inputs/meta_analysis_v1/HFD_primary105_meta_effects_primaryQC_v1.tsv.gz'
    for f in [core_path, loo_path, effects_path]:
        if not f.is_file():
            parser.error(f'Missing required input file: {f}')

    core_rows = list(read_tsv(core_path))
    core = {r['Geneid']: r for r in core_rows}
    if len(core) != 59 or len(core_rows) != 59:
        raise ValueError(f'Expected 59 distinct original core genes, got {len(core)} IDs/{len(core_rows)} rows')

    # Reproduce GSE-level inverse-variance summary using the exact primary QC inputs.
    # Streaming means memory scales with 59 candidate genes, not all >1.5M effects.
    ivw = defaultdict(lambda: [0.0, 0.0, 0])
    candidate_rows = 0
    for r in read_tsv(effects_path):
        gene = r['Geneid']
        if gene not in core:
            continue
        y = number(r, 'log2FoldChange')
        v = number(r, 'sampling_variance')
        if v <= 0:
            raise ValueError(f'Nonpositive variance for {gene}')
        item = ivw[(gene, r['GSE'])]
        item[0] += y / v
        item[1] += 1.0 / v
        item[2] += 1
        candidate_rows += 1

    summaries = defaultdict(dict)
    contrasts = defaultdict(dict)
    for (gene, gse), (yw, w, n) in ivw.items():
        summaries[gene][gse] = yw / w
        contrasts[gene][gse] = n

    for gene, r in core.items():
        if gene not in summaries:
            raise ValueError(f'No meta-analysis input effects for core gene {gene}')
        actual_n = len(summaries[gene])
        expected_n = int(float(r['N_independent_GSE']))
        if actual_n != expected_n:
            raise ValueError(f'{gene}: GSE count differs (input={actual_n}; core={expected_n})')
        orig_cons = sum(sign(e) == sign(number(r, 'pooled_log2FC'))
                        for e in summaries[gene].values()) / actual_n
        stored_cons = number(r, 'GSE_balanced_direction_consistency')
        if not math.isclose(orig_cons, stored_cons, abs_tol=1e-11):
            raise ValueError(f'{gene}: original directional consistency not reproduced '
                             f'({orig_cons} vs {stored_cons})')

    detailed = []
    seen = set()
    for r in read_tsv(loo_path):
        gene = r['Geneid']
        if gene not in core:
            raise ValueError(f'Unexpected LOO gene in detail: {gene}')
        dropped = r['Dropped_GSE']
        key = (gene, dropped)
        if key in seen:
            raise ValueError(f'Duplicate LOO gene/GSE pair: {key}')
        seen.add(key)
        if dropped not in summaries[gene]:
            raise ValueError(f'{gene}: dropped GSE {dropped} is absent from input')
        remaining = {g: e for g, e in summaries[gene].items() if g != dropped}
        n_remain = len(remaining)
        if n_remain != int(float(r['Remaining_GSE'])):
            raise ValueError(f'{gene}/{dropped}: inconsistent remaining GSE count')
        n_contrasts = sum(contrasts[gene][g] for g in remaining)
        if n_contrasts != int(float(r['Remaining_contrasts'])):
            raise ValueError(f'{gene}/{dropped}: inconsistent remaining contrast count')
        status = r['LOO_status'].upper()
        out = {
            'Geneid': gene,
            'Dropped_GSE': dropped,
            'N_full_GSE': len(summaries[gene]),
            'N_remaining_GSE': n_remain,
            'Meets_original_min_10_GSE': int(n_remain >= 10),
            'LOO_status': status,
            'LOO_pooled_log2FC': '',
            'LOO_p_value_nominal': '',
            'LOO_PI_lower_95': '',
            'LOO_PI_upper_95': '',
            'LOO_PI_excludes_zero': '',
            'LOO_GSE_balanced_consistency': '',
            'LOO_consistency_ge_0_80': '',
            'LOO_PI_and_consistency_both': '',
            'LOO_nominal_p_lt_0_05': '',
        }
        if status == 'PASS':
            mu = number(r, 'LOO_pooled_log2FC')
            se = number(r, 'LOO_SE')
            tau = number(r, 'LOO_tau2_total')
            p = number(r, 'LOO_p_value')
            if se < 0 or tau < -1e-10:
                raise ValueError(f'{gene}/{dropped}: negative SE or variance')
            sd = math.sqrt(se * se + max(0.0, tau))
            low, high = mu - Z95 * sd, mu + Z95 * sd
            pi_ok = low > 0 or high < 0
            consistency = sum(sign(e) == sign(mu) for e in remaining.values()) / n_remain
            cons_ok = consistency >= 0.8 - 1e-12
            out.update({
                'LOO_pooled_log2FC': f'{mu:.12g}',
                'LOO_p_value_nominal': f'{p:.12g}',
                'LOO_PI_lower_95': f'{low:.12g}',
                'LOO_PI_upper_95': f'{high:.12g}',
                'LOO_PI_excludes_zero': int(pi_ok),
                'LOO_GSE_balanced_consistency': f'{consistency:.12g}',
                'LOO_consistency_ge_0_80': int(cons_ok),
                'LOO_PI_and_consistency_both': int(pi_ok and cons_ok),
                'LOO_nominal_p_lt_0_05': int(p < 0.05),
            })
            full_sign = sign(number(core[gene], 'pooled_log2FC'))
            if sign(mu) != full_sign:
                raise ValueError(f'{gene}/{dropped}: pooled LOO direction flip, inconsistent with reported all-stable result')
        elif status != 'FAIL':
            raise ValueError(f'Unexpected LOO status: {status}')
        detailed.append(out)

    by_gene = []
    for gene in sorted(core):
        vals = [r for r in detailed if r['Geneid'] == gene]
        expected = len(summaries[gene])
        if len(vals) != expected:
            raise ValueError(f'{gene}: expected {expected} LOO rows, got {len(vals)}')
        passed = [r for r in vals if r['LOO_status'] == 'PASS']
        if not passed:
            raise ValueError(f'{gene}: zero successful LOO refits')
        n_pi = sum(r['LOO_PI_excludes_zero'] == 1 for r in passed)
        n_cons = sum(r['LOO_consistency_ge_0_80'] == 1 for r in passed)
        n_both = sum(r['LOO_PI_and_consistency_both'] == 1 for r in passed)
        n_nominal = sum(r['LOO_nominal_p_lt_0_05'] == 1 for r in passed)
        n_with10 = sum(r['Meets_original_min_10_GSE'] == 1 for r in vals)
        by_gene.append({
            'Geneid': gene,
            'N_full_GSE': expected,
            'N_LOO_total': len(vals),
            'N_LOO_fit_PASS': len(passed),
            'N_LOO_PI_excludes_zero': n_pi,
            'N_LOO_consistency_ge_0_80': n_cons,
            'N_LOO_PI_and_consistency_both': n_both,
            'N_LOO_nominal_p_lt_0_05': n_nominal,
            'N_LOO_remaining_GSE_ge_10': n_with10,
            'All_LOO_fit_PASS': int(len(passed) == expected),
            'All_LOO_PI_pass': int(n_pi == expected),
            'All_LOO_consistency_pass': int(n_cons == expected),
            'All_LOO_PI_and_consistency_pass': int(n_both == expected),
            'All_LOO_nominal_p_lt_0_05': int(n_nominal == expected),
            'All_LOO_keep_10_GSE_eligibility': int(n_with10 == expected),
            'Min_LOO_directional_consistency': f'{min(float(r["LOO_GSE_balanced_consistency"]) for r in passed):.8g}',
        })

    n_loo = len(detailed)
    n_fail = sum(r['LOO_status'] == 'FAIL' for r in detailed)
    n_coverage_drop = sum(r['Meets_original_min_10_GSE'] == 0 for r in detailed)
    result = [
        'IJMS R1 Major Comment 1: leave-one-contributing-GSE-out PI/consistency audit',
        f'Original core genes={len(core)}',
        f'Original gene-contrast rows in primary QC input for these genes={candidate_rows}',
        f'Original existing LOO refits={n_loo}',
        f'Existing LOO refit failures={n_fail}',
        f'LOO refits with remaining gene coverage below original >=10 GSE eligibility={n_coverage_drop}',
        f'Genes keeping 95% PI exclusion in every LOO refit={sum(r["All_LOO_PI_pass"] for r in by_gene)}',
        f'Genes keeping consistency>=0.80 in every LOO refit={sum(r["All_LOO_consistency_pass"] for r in by_gene)}',
        f'Genes keeping BOTH PI exclusion and consistency>=0.80 in every LOO refit={sum(r["All_LOO_PI_and_consistency_pass"] for r in by_gene)}',
        f'Genes with nominal p<0.05 in every LOO refit={sum(r["All_LOO_nominal_p_lt_0_05"] for r in by_gene)}',
        '',
        'WARNING: These are NOT complete three-criterion core classifications. The original pooled BH-FDR ',
        'cannot be reevaluated from LOO fits for only the 59 candidates; genome-wide LOO refits would be needed.',
        'LOO fits for genes with exactly 10 contributing GSEs retain only 9 GSEs and would fail original ',
        'meta-analysis eligibility. Such fits are exploratory influence diagnostics, not new primary core assignments.',
        'The analysis reuses stored LOO multilevel fits; it does not perform new model fits.',
    ]

    outdir = (args.outdir or root / '00_metadata/analysis_outputs/revision_R1C1_LOO_criteria').resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(outdir / 'loo_pi_consistency_detail_v1.tsv', detailed, list(detailed[0]))
    write_tsv(outdir / 'loo_pi_consistency_by_gene_v1.tsv', by_gene, list(by_gene[0]))
    (outdir / 'audit_summary_v1.txt').write_text('\n'.join(result) + '\n', encoding='utf-8')
    print('\n'.join(result))
    print(f'\nOutputs: {outdir}')


if __name__ == '__main__':
    main()
