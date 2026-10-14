#!/usr/bin/env python3
"""Cross-check final frozen moderator coverage against Batch06/Batch11 audit conventions.

READ ONLY: never alters analysis inputs, tables, or models.
Compares by Contrast_id and distinguishes *descriptive nonmissing* coverage
from Batch11 'known' and 'inclusive-ready' coverage.
"""
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

FIELDS = ('hfd_duration_weeks', 'hfd_fat_kcal_percent', 'sex', 'strain', 'control_diet_type')
MISSING = {'', 'na', 'n/a', 'nan', 'none', 'null', 'unknown', 'not_reported',
           'not reported', 'unresolved', 'pending', 'missing', '?', 'not available',
           'not determined', 'not specified', 'unspecified'}
STRAINS = {'C57BL/6J', 'C57BL/6', 'C57BL/6NHsd', 'C57BL/6N',
           'C57BL/6NTac', '129S1/SvImJ', 'FVB.129P2'}
CONTROLS = {'CHOW', 'LFD', 'NORMAL_DIET', 'CONTROL_DIET_UNSPECIFIED'}


def read_table(path):
    with path.open(newline='', encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def descriptive_present(value, var):
    v = str(value or '').strip()
    if v.lower() in MISSING:
        return False
    if var in ('hfd_duration_weeks', 'hfd_fat_kcal_percent'):
        try:
            return math.isfinite(float(v))
        except ValueError:
            return False
    return True


def final_effective(row, var):
    """Reproduce script 65 Batch11 value origin rules."""
    if row.get(f'{var}_adjudication_status', '').strip() != 'NOT_REVIEWED':
        return row.get(f'{var}_final', '').strip(), 'MANUAL'
    if row.get(f'{var}_geo_auto_status', '').strip() == 'RESOLVED':
        return row.get(f'{var}_geo_auto_value', '').strip(), 'GEO_AUTO'
    return '', 'NONE'


def inclusive_ready(var, val):
    if not val or val in {'UNKNOWN', 'UNRESOLVED', 'NOT_REPORTED', 'NA', 'N/A'}:
        return False
    if var == 'hfd_duration_weeks':
        try:
            return float(val) > 0
        except ValueError:
            return False
    if var == 'hfd_fat_kcal_percent':
        try:
            return 0 < float(val) <= 100
        except ValueError:
            return False
    if var == 'sex':
        return val in {'MALE', 'FEMALE'}
    if var == 'strain':
        return val in STRAINS
    if var == 'control_diet_type':
        return val in CONTROLS
    return False


def as_int(value):
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return None


def summarize(report, label):
    if not report.is_file():
        return {}, f'{label}: not present ({report})'
    rr = {r['Moderator']: r for r in read_table(report)}
    return rr, f'{label}: {report}'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path.cwd())
    a = p.parse_args()
    root = a.root.resolve()
    frozen = root/'00_metadata/analysis_inputs/frozen_v1/HFD_primary105_contrast_analysis_input_audit_v1_1.tsv'
    master = root/'00_metadata/contrasts/moderators/manual_adjudication/HFD_primary105_moderator_manual_adjudication_batch11_v1.tsv'
    covdir = root/'00_metadata/contrasts/moderators/manual_adjudication'
    cov06, status06 = summarize(covdir/'HFD_primary105_moderator_coverage_after_batch06_v1_1.tsv', 'Batch06 report')
    cov11, status11 = summarize(covdir/'HFD_primary105_moderator_coverage_after_batch11_v1.tsv', 'Batch11 report')
    if not frozen.exists() or not master.exists():
        raise SystemExit(f'Missing frozen or Batch11 master:\n{frozen}\n{master}')
    frows = read_table(frozen)
    mrows = read_table(master)
    f = {r['Contrast_id'].strip(): r for r in frows}
    m = {r['Contrast_id'].strip(): r for r in mrows}
    if len(frows) != 105 or len(mrows) != 105 or len(f) != 105 or len(m) != 105 or set(f) != set(m):
        raise SystemExit('Expected exactly 105 matching, unique Contrast_id in both inputs.')
    if len({r['GSE'] for r in frows}) != 64:
        raise SystemExit('Expected 64 GSEs in frozen table.')

    outdir = root/'00_metadata/analysis_outputs/revision_R1C6_moderator_confounds'
    outdir.mkdir(parents=True, exist_ok=True)
    cols = ['Moderator','Frozen_descriptive_nonmissing','Batch11_final_descriptive_nonmissing',
            'Batch11_effective_inclusive_ready','Batch06_report_known','Batch06_report_inclusive_ready',
            'Batch11_report_known','Batch11_report_inclusive_ready','Frozen_alias_vs_final_mismatch']
    summary = []
    diffs = []
    for var in FIELDS:
        present = final_present = inc = alias_diff = 0
        for cid in sorted(f):
            x = f[cid]; r = m[cid]
            if x['GSE'].strip() != r['GSE'].strip():
                raise SystemExit(f'{cid}: GSE mismatch')
            frozen_val = x.get(var, '').strip()
            final_val = r.get(f'{var}_final', '').strip()
            effect, source = final_effective(r,var)
            a_present = descriptive_present(frozen_val, var)
            b_present = descriptive_present(final_val, var)
            c_ready = inclusive_ready(var,effect)
            alias_mismatch = frozen_val != final_val
            present += a_present; final_present += b_present
            inc += c_ready; alias_diff += alias_mismatch
            if a_present != c_ready or alias_mismatch:
                diffs.append({
                    'Moderator':var, 'Contrast_id':cid,'GSE':x['GSE'],
                    'Frozen_alias_value': frozen_val,
                    'Batch11_final_value': final_val,
                    'Batch11_effective_value': effect,
                    'Effective_source':source,
                    'Adjudication_status':r.get(f'{var}_adjudication_status',''),
                    'GEO_auto_status':r.get(f'{var}_geo_auto_status',''),
                    'Frozen_descriptive_nonmissing':int(a_present),
                    'Batch11_inclusive_ready':int(c_ready),
                    'Alias_vs_final_differ':int(alias_mismatch),
                })
        summary.append({
            'Moderator':var,
            'Frozen_descriptive_nonmissing':present,
            'Batch11_final_descriptive_nonmissing':final_present,
            'Batch11_effective_inclusive_ready':inc,
            'Batch06_report_known':cov06.get(var,{}).get('Known_contrasts','NA'),
            'Batch06_report_inclusive_ready':cov06.get(var,{}).get('Inclusive_ready_contrasts','NA'),
            'Batch11_report_known':cov11.get(var,{}).get('Known_contrasts','NA'),
            'Batch11_report_inclusive_ready':cov11.get(var,{}).get('Inclusive_ready_contrasts','NA'),
            'Frozen_alias_vs_final_mismatch':alias_diff,
        })
    output1 = outdir/'coverage_crosscheck_summary_v1.tsv'
    output2 = outdir/'coverage_crosscheck_contrast_differences_v1.tsv'
    with output1.open('w', newline='', encoding='utf-8') as fh:
        w=csv.DictWriter(fh, fieldnames=cols,delimiter='\t'); w.writeheader(); w.writerows(summary)
    diffcols = ['Moderator','Contrast_id','GSE','Frozen_alias_value','Batch11_final_value',
                'Batch11_effective_value','Effective_source','Adjudication_status','GEO_auto_status',
                'Frozen_descriptive_nonmissing','Batch11_inclusive_ready','Alias_vs_final_differ']
    with output2.open('w', newline='', encoding='utf-8') as fh:
        w=csv.DictWriter(fh, fieldnames=diffcols,delimiter='\t');w.writeheader();w.writerows(diffs)
    print('===== R1C6 moderator coverage cross-check =====')
    print('Input frozen:',frozen)
    print('Input master:',master)
    print(status06)
    print(status11)
    print('\t'.join(cols))
    for r in summary:
        print('\t'.join(str(r[k]) for k in cols))
    print(f'\nFound {len(diffs)} moderator × contrast differences between frozen descriptive presence and Batch11 inclusive readiness (or raw alias mismatch).')
    print('Summary:',output1)
    print('Per-contrast differences:',output2)
    print('Interpretation: Frozen descriptive/nonmissing and Batch11 inclusive-ready are DIFFERENT denominators.')
    print('Do not overwrite frozen analyses based on differences in these denominators.')

if __name__=='__main__':
    main()
