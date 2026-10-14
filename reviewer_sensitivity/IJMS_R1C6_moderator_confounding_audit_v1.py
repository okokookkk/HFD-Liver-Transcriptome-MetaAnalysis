#!/usr/bin/env python3
"""Reviewer 1, major comment 6: independent-GSE-aware descriptive moderator audit.

Input: frozen, publication-level 105-contrast design table (same source as original Table S2).
All associations below are descriptive, not gene-level regression or causal estimates.
Only study-level summaries are used for pairwise associations, to prevent a GSE
with multiple contrasts from being counted as several independent experiments.
"""
import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

REL = '00_metadata/analysis_inputs/frozen_v1/HFD_primary105_contrast_analysis_input_audit_v1_1.tsv'
OUTREL = '00_metadata/analysis_outputs/revision_R1C6_moderator_confounds'
NUMERIC = ('hfd_duration_weeks', 'hfd_fat_kcal_percent')
CATEGORICAL = ('sex', 'strain', 'control_diet_type')
FIELDS = NUMERIC + CATEGORICAL
MISSING = {'', 'na', 'n/a', 'nan', 'none', 'null', 'unknown', 'not_reported',
           'not reported', 'unresolved', 'pending', 'missing', '?', 'not available',
           'not determined', 'not specified', 'unspecified'}


def clean(value, name):
    v = str(value if value is not None else '').strip()
    if v.lower() in MISSING:
        return None
    if name in NUMERIC:
        try:
            x = float(v)
            return x if math.isfinite(x) else None
        except ValueError:
            return None
    v = ' '.join(v.split()).upper()
    if name == 'sex':
        return {'M': 'MALE', 'F': 'FEMALE', 'MALES': 'MALE',
                'FEMALES': 'FEMALE', 'BOTH': 'BOTH/MIXED',
                'MIXED': 'BOTH/MIXED', 'M/F': 'BOTH/MIXED'}.get(v, v)
    return v


def qtile(xs, prob):
    xs = sorted(xs)
    if not xs:
        return float('nan')
    r = (len(xs) - 1) * prob
    l = int(math.floor(r))
    u = int(math.ceil(r))
    return xs[l] * (u-r) + xs[u] * (r-l) if l != u else xs[l]


def rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ret = [0.] * len(xs)
    j = 0
    while j < len(order):
        k = j + 1
        while k < len(order) and xs[order[k]] == xs[order[j]]:
            k += 1
        midrank = (j + 1 + k) / 2.
        for z in order[j:k]:
            ret[z] = midrank
        j = k
    return ret


def corr(x, y):
    n = len(x)
    if n < 3:
        return float('nan')
    mx = sum(x)/n; my = sum(y)/n
    dx = [v-mx for v in x]; dy = [v-my for v in y]
    den = math.sqrt(sum(v*v for v in dx)*sum(v*v for v in dy))
    return sum(a*b for a, b in zip(dx, dy))/den if den > 0 else float('nan')


def spearman(x, y):
    return corr(rank(x), rank(y))


def cramer_v(pairs):
    if len(pairs) < 3:
        return float('nan')
    a = Counter(x for x, _ in pairs)
    b = Counter(y for _, y in pairs)
    if len(a) < 2 or len(b) < 2:
        return float('nan')
    both = Counter(pairs)
    n = len(pairs)
    chisq = sum((both.get((x,y),0)-a[x]*b[y]/n)**2/(a[x]*b[y]/n)
                for x in a for y in b)
    # Raw Cramer's V, descriptive only; small/sparse groups warrant caution.
    return math.sqrt(chisq/(n*min(len(a)-1, len(b)-1)))


def write_tsv(path, cols, rows):
    with path.open('w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow({k: f'{v:.6g}' if isinstance(v,float) and math.isfinite(v)
                        else ('NA' if isinstance(v,float) else v) for k,v in row.items()})


def paired_studies(studies, a, b):
    """For each GSE, keep only contrasts where BOTH covariates are observed."""
    out=[]
    for gse, rs in studies.items():
        pairs = [(r[a],r[b]) for r in rs if r[a] is not None and r[b] is not None]
        if pairs:
            out.append((gse, pairs))
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path.cwd())
    args=parser.parse_args()
    root=args.root.expanduser().resolve()
    source=root/REL
    if not source.is_file():
        raise SystemExit('Required frozen metadata file not found:\n  '+str(source))
    with source.open(encoding='utf-8-sig', newline='') as fh:
        reader=csv.DictReader(fh,delimiter='\t')
        required=['Contrast_id','GSE',*FIELDS]
        missing=[k for k in required if k not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit('Missing required columns: '+', '.join(missing))
        rows=[]
        for r in reader:
            rows.append({'Contrast_id':r['Contrast_id'].strip(),
                         'GSE':r['GSE'].strip(),
                         **{k:clean(r.get(k),k) for k in FIELDS}})
    if len(rows)!=105 or len(set(r['Contrast_id'] for r in rows))!=105:
        raise SystemExit(f'Expected 105 unique contrasts; found {len(rows)} rows')
    studies=defaultdict(list)
    for r in rows:
        studies[r['GSE']].append(r)
    if len(studies)!=64:
        raise SystemExit(f'Expected 64 GSEs; found {len(studies)}')

    outdir=root/OUTREL
    outdir.mkdir(parents=True,exist_ok=True)
    coverage=[]; overview=[]
    for k in FIELDS:
        nc=sum(r[k] is not None for r in rows)
        ng=sum(any(r[k] is not None for r in rr) for rr in studies.values())
        multiple=sum(len({r[k] for r in rr if r[k] is not None})>1 for rr in studies.values())
        coverage.append(dict(Moderator=k,N_contrasts_known=nc,Total_contrasts=105,
                             N_GSE_with_known=ng,Total_GSE=64,
                             N_GSE_with_multiple_recorded_values=multiple))
    write_tsv(outdir/'metadata_coverage_v1.tsv', list(coverage[0]),coverage)

    for gse, rr in sorted(studies.items()):
        d={'GSE':gse,'N_contrasts':len(rr)}
        for k in FIELDS:
            vals=sorted({r[k] for r in rr if r[k] is not None})
            d[k+'_n_unique']=len(vals)
            d[k+'_observed_values']=' | '.join(map(str,vals))
        overview.append(d)
    write_tsv(outdir/'gse_covariate_inventory_v1.tsv',list(overview[0]),overview)

    nn=[]
    for i,a in enumerate(NUMERIC):
        for b in NUMERIC[i+1:]:
            ps=paired_studies(studies,a,b)
            # Only pair-complete contrasts per GSE; each GSE contributes one median pair.
            x=[median([p[0] for p in vv]) for _,vv in ps]
            y=[median([p[1] for p in vv]) for _,vv in ps]
            xc=[p[0] for _,vv in ps for p in vv]
            yc=[p[1] for _,vv in ps for p in vv]
            nn.append(dict(Moderator_A=a,Moderator_B=b,
                  N_pair_complete_contrasts=len(xc),N_independent_GSE=len(ps),
                  Spearman_rho_GSE_equal_weight=spearman(x,y),
                  Spearman_rho_contrast_level_descriptive_only=spearman(xc,yc),
                  Note='GSE-level medians use same pair-complete contrasts for both variables'))
    write_tsv(outdir/'duration_fat_spearman_v1.tsv',list(nn[0]),nn)

    groups=[]; diagnostics=[]
    for a in NUMERIC:
        for b in CATEGORICAL:
            ps=paired_studies(studies,a,b)
            gs=defaultdict(list)
            excluded=[]
            for gse,vv in ps:
                cats={p[1] for p in vv}
                if len(cats)!=1:
                    excluded.append(gse)
                    continue
                gs[next(iter(cats))].append(median([p[0] for p in vv]))
            diagnostics.append(dict(Continuous_moderator=a,Categorical_moderator=b,
                N_GSE_pair_complete=len(ps),N_GSE_excluded_multiple_categories=len(excluded),
                N_GSE_retained=sum(map(len,gs.values())),N_categories=len(gs)))
            for cat,xx in sorted(gs.items()):
                groups.append(dict(Continuous_moderator=a,Categorical_moderator=b,
                    Category=cat,N_independent_GSE=len(xx),Median=median(xx),
                    Q1=qtile(xx,.25),Q3=qtile(xx,.75),
                    Interpretation_note=('SMALL_GROUP' if len(xx)<3 else 'DESCRIPTIVE_ONLY')))
    write_tsv(outdir/'numeric_by_category_v1.tsv',list(groups[0]),groups)
    write_tsv(outdir/'numeric_by_category_diagnostics_v1.tsv',list(diagnostics[0]),diagnostics)

    cc=[]
    for i,a in enumerate(CATEGORICAL):
        for b in CATEGORICAL[i+1:]:
            ps=paired_studies(studies,a,b)
            vals=[]; skipped=0
            for gse,vv in ps:
                aa={p[0] for p in vv}; bb={p[1] for p in vv}
                if len(aa)!=1 or len(bb)!=1:
                    skipped+=1
                    continue
                vals.append((next(iter(aa)),next(iter(bb))))
            cc.append(dict(Categorical_A=a,Categorical_B=b,
                N_GSE_pair_complete=len(ps),N_GSE_retained=len(vals),
                N_GSE_excluded_multiple_categories=skipped,
                N_categories_A=len(set(a for a,b in vals)),
                N_categories_B=len(set(b for a,b in vals)),
                Cramers_V_descriptive=cramer_v(vals),
                Note='Sparse categories and within-GSE mixed designs limit interpretation'))
    write_tsv(outdir/'categorical_pair_cramers_v_v1.tsv',list(cc[0]),cc)

    print('===== Reviewer 1 Comment 6: Moderator confounding diagnostic =====')
    print('Source:',source)
    print(f'GSEs={len(studies)}; contrasts={len(rows)}')
    print('\nCovariate coverage:')
    for r in coverage:
        print(f"  {r['Moderator']}: {r['N_contrasts_known']}/105 contrasts; "
              f"{r['N_GSE_with_known']}/64 GSEs; "
              f"within-GSE multiple values in {r['N_GSE_with_multiple_recorded_values']} GSEs")
    print('\nDuration vs fat percentage (GSE-balanced):')
    for r in nn:
        print(f"  Pair-complete: {r['N_pair_complete_contrasts']} contrasts / "
              f"{r['N_independent_GSE']} GSEs; Spearman rho={r['Spearman_rho_GSE_equal_weight']:.3f}")
    print('\nDuration by sex, strain, and control diet:')
    for d in diagnostics:
        if d['Continuous_moderator']=='hfd_duration_weeks':
            print(f"  {d['Categorical_moderator']}: {d['N_GSE_retained']} GSEs, "
                  f"{d['N_categories']} categories; {d['N_GSE_excluded_multiple_categories']} mixed-GSE excluded")
    print('\nOutputs:',outdir)
    print('IMPORTANT: These are DESCRIPTIVE cross-study associations only. Unknown age,')
    print('platform, study design, and metabolic phenotype are NOT controlled by this audit.')
    print('Multiple GSE contrasts are not treated as independent; each GSE contributes')
    print('at most one summary to an association estimate. No causal conclusions are possible.')

if __name__=='__main__':
    main()
