#!/usr/bin/env python3
"""Reviewer 1 Major Comment 5: descriptive heterogeneity component summary.

Reads the frozen primary gene-level multilevel REML meta-analysis output.
The 59-core assignment is NOT refitted or inferred from heterogeneity alone.
"""
import argparse
from pathlib import Path
import pandas as pd

REL = '05_meta/gene_level_v1_2/HFD_gene_multilevel_meta_primary_v1_2.tsv.gz'
COLS = ['Geneid', 'pooled_FDR_BH', 'stable_core_preLOO',
        'tau2_GSE', 'tau2_contrast_within_GSE', 'tau2_total',
        'I2_GSE_percent', 'I2_contrast_within_GSE_percent', 'I2_total_percent']
METRICS = COLS[3:]

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path.cwd())
    args = p.parse_args()
    root = args.root.expanduser().resolve()
    source = root / REL
    if not source.is_file():
        raise SystemExit(f'Missing required source: {source}')
    df = pd.read_csv(source, sep='\t', compression='gzip', usecols=COLS)
    if df.Geneid.duplicated().any():
        raise SystemExit('ERROR: duplicate Geneid rows in primary gene-level output')
    for c in ['pooled_FDR_BH'] + METRICS:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    if df[METRICS].isna().any().any():
        raise SystemExit('ERROR: missing/non-numeric values in heterogeneity metrics')
    sig = df.pooled_FDR_BH < 0.05
    core = df.stable_core_preLOO.astype(str).str.upper() == 'YES'
    if (core & ~sig).any():
        raise SystemExit('ERROR: core contains FDR-nonsignificant genes')
    groups = {
        'All meta-analyzed genes': df,
        'Pooled FDR <0.05 (including core)': df[sig],
        'Pooled FDR <0.05, excluding core': df[sig & ~core],
        'Original 59-gene core': df[core],
    }
    rows = []
    for name, sub in groups.items():
        for metric in METRICS:
            x = sub[metric]
            rows.append({
                'Gene_group': name,
                'N_genes': len(sub),
                'Measure': metric,
                'Median': x.median(),
                'Q1': x.quantile(0.25),
                'Q3': x.quantile(0.75),
                'Proportion_zero_percent': 100 * (x == 0).mean(),
            })
    outdir = root / '00_metadata/analysis_outputs/revision_R1C5_heterogeneity'
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / 'heterogeneity_component_summary_v1.tsv'
    pd.DataFrame(rows).to_csv(out, sep='\t', index=False, float_format='%.6g')

    print('===== R1 Major Comment 5: Heterogeneity components =====')
    print('Original gene-level results:', source)
    print(f'Genes total={len(df)}; pooled FDR<0.05={sig.sum()}; original core={core.sum()}')
    if (len(df), int(sig.sum()), int(core.sum())) != (18234, 6920, 59):
        print('WARNING: primary analysis counts differ from manuscript; verify source version.')
    for name, sub in groups.items():
        print(f'\n{name} (n={len(sub)})')
        for metric in METRICS:
            x = sub[metric]
            print(f'  {metric}: median={x.median():.3f}, IQR={x.quantile(.25):.3f}–{x.quantile(.75):.3f}; 0-estimate={100*(x==0).mean():.1f}%')
    print(f'\nSaved: {out}')
    print('NOTE: lower heterogeneity in a PI-selected core is partly selection-induced and not an independent validation.')

if __name__ == '__main__':
    main()
