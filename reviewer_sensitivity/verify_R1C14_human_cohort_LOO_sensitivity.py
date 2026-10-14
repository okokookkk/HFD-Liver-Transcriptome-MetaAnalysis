#!/usr/bin/env python3
"""Exploratory Reviewer 1 Comment 14 human-cohort-exclusion sensitivity.
Inputs: submitted SupplementaryTables_v20260922.xlsx, sheets S14,S16.
Uses fixed-effect inverse-variance pooling for 2 remaining cohorts, with
mouse-human same-sign endpoint regressed on +0.10 GSE directional consistency.
The FE LOO effect is NOT a 2-cohort REML re-fit; do not label as such.
"""
from zipfile import ZipFile
from lxml import etree
import numpy as np
import statsmodels.api as sm
import csv
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--workbook', type=Path, required=True,
                    help='Exact submitted SupplementaryTables_v20260922.xlsx (S14 and S16)')
parser.add_argument('--output', type=Path, default=Path('outputs/HFD_R1C14_human_cohort_exclusion_sensitivity.csv'))
args = parser.parse_args()
src = str(args.workbook.expanduser().resolve())
if not Path(src).is_file():
    parser.error(f'Missing workbook: {src}')
NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

def get_sheet_rows(z,sheet_idx,cols):
    # Parse source XML in place; ignore formatting and unneeded columns.
    contents=z.read(f'xl/worksheets/sheet{sheet_idx}.xml')
    wanted=set(cols)
    from io import BytesIO
    for ev,node in etree.iterparse(BytesIO(contents),events=('end',),tag=NS+'row'):
        vals={}
        for cell in node:
            addr=cell.get('r','')
            col=''.join(x for x in addr if x.isalpha())
            if col not in wanted: continue
            value=cell.find(NS+'v')
            if value is None: continue
            vals[col]=value.text
        yield vals
        node.clear()

with ZipFile(src) as z:
    s14=get_sheet_rows(z,15,['A','B','C','F','G','J','K'])
    h=next(s14);print('Human cohort columns:',h)
    effects={}
    for r in s14:
        try:
            key=str(int(float(r['A'])))
            betas=np.array([float(r[c]) for c in ('B','F','J')])
            ses=np.array([float(r[c]) for c in ('C','G','K')])
            if np.isfinite(betas).all() and np.isfinite(ses).all() and (ses>0).all():
                effects[key]=(betas,ses)
        except (ValueError,KeyError): pass
    s16=get_sheet_rows(z,17,['A','F','H','K','M'])
    hm=next(s16);print('mouse-human cols:',hm)
    ortho=[]
    for r in s16:
        try:
            key=str(int(float(r['A'])))
            if key not in effects:continue
            mouse_beta=float(r['F']); cons=float(r['H']); pooled=float(r['K'])
            if 0<=cons<=1 and np.isfinite(mouse_beta) and np.isfinite(pooled):
                ortho.append((key,mouse_beta,cons,pooled))
        except (ValueError,KeyError): pass

keys=[k for k,_,_,_ in ortho];mouse_beta=np.array([row[1] for row in ortho]);cons=np.array([row[2] for row in ortho]);hre=np.array([row[3] for row in ortho]);
betas=np.stack([effects[k][0] for k in keys]);ses=np.stack([effects[k][1] for k in keys]);

print('n='+str(len(keys)))

def glm_result(name,beta_h,valid=None):
    if valid is None:valid=np.ones(len(keys),bool)
    valid=valid & np.isfinite(beta_h) & (beta_h!=0)&(mouse_beta!=0)
    x=cons[valid]/.1
    y=(np.sign(beta_h[valid])==np.sign(mouse_beta[valid])).astype(int)
    res=sm.GLM(y,sm.add_constant(x),family=sm.families.Binomial()).fit()
    c=float(res.params[1]);se=float(res.bse[1]);o=np.exp(c);ci=np.exp(np.array([c-1.96*se,c+1.96*se]));p=float(res.pvalues[1]);
    return [name,int(valid.sum()),int(y.sum()),round(float(y.mean()),6),float(o),float(ci[0]),float(ci[1]),p]

w=1/ses**2
allfe=(w*betas).sum(1)/w.sum(1)
res=[glm_result('Primary 3-cohort RE',hre),glm_result('3-cohort FE',allfe)]
for j,name in enumerate(['GSE126848','GSE130970','GSE135251']):
    keep=[i for i in range(3) if i!=j]
    bw=betas[:,keep];ww=w[:,keep]
    pooled=(ww*bw).sum(1)/ww.sum(1)
    res.append(glm_result('Exclude '+name+' (2-cohort FE)',pooled))

out=str(args.output.expanduser().resolve())
Path(out).parent.mkdir(parents=True, exist_ok=True)
with open(out,'w',newline='') as f:
    cw=csv.writer(f);cw.writerow(['Comparison','N_genes','N_same_direction','Same_direction_fraction','OR_per_10pp','CI95_lower','CI95_upper','P_value']);cw.writerows(res)
for r in res: print(r)
print('Created',out)
