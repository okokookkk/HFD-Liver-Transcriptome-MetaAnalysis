#!/usr/bin/env python3
"""Inventory crucial frozen metadata and selected upstream analytical data; does not run models."""
import argparse
from pathlib import Path

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, default=Path('.'))
parser.add_argument('--strict', action='store_true', help='Nonzero exit status if required scientific inputs missing')
args=parser.parse_args()
root=args.root.expanduser().resolve()
required=[
'00_metadata/frozen_search/HFD_primary_analysis_GSE_frozen_v1.txt',
'00_metadata/analysis_inputs/frozen_v1/HFD_primary105_contrast_analysis_input_audit_v1_1.tsv',
'00_metadata/technical_harmonization/frozen_v1/HFD_primary_ENA_R1R2_download_manifest_frozen_v1.tsv',
'00_metadata/contrasts/moderators/manual_adjudication/HFD_primary105_moderator_manual_adjudication_batch11_v1.tsv',
]
upstream=[
'00_metadata/analysis_inputs/meta_analysis_v1/HFD_primary105_meta_effects_primaryQC_v1.tsv.gz',
'05_meta/gene_level_v1_2/HFD_gene_multilevel_meta_primary_v1_2.tsv.gz',
'05_meta/stable_core_v1/HFD_stable_core_final_v1.tsv',
'07_human_validation/gene_meta_v1/human_MASLD_three_cohort_random_effects_meta_v1.tsv.gz',
]
print('Project root:',root)
miss=[]
for p in required+upstream:
    ok=(root/p).is_file() and (root/p).stat().st_size>0
    print(('OK     ' if ok else 'MISSING'),p)
    if not ok:miss.append(p)
print('Missing:',len(miss),'of',len(required+upstream))
if args.strict and miss:raise SystemExit(2)
