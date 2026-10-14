#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

ACQ = Path("00_metadata/data_acquisition")
OUTDIR = Path("00_metadata/sample_screening")
OUTDIR.mkdir(parents=True, exist_ok=True)

INFILE = ACQ/"Data_Acquisition_Master_v1.tsv"
OUT = OUTDIR/"primary107_sample_review_pack_v1.txt"
SUMMARY = OUTDIR/"primary107_sample_review_summary_v1.tsv"

df = pd.read_csv(INFILE, sep="\t", dtype=str).fillna("")

display_cols = [c for c in [
    "GSM","Sample_title","Source_name","Tissue","Diet","Genotype","Sex","Age","Strain",
    "Treatment","Library_strategy","Raw_availability","N_SRR","Sample_review_status_v1"
] if c in df.columns]

with OUT.open("w") as f:
    for idx, (gse, sub) in enumerate(df.groupby("GSE", sort=True), start=1):
        f.write("#"*100+"\n")
        f.write(f"### {gse}   [{idx}/107]\n")
        f.write("#"*100+"\n")
        counts = sub["Sample_review_status_v1"].value_counts().to_dict()
        f.write("Status counts: "+" | ".join(f"{k}={v}" for k,v in sorted(counts.items()))+"\n\n")
        f.write(sub[display_cols].to_string(index=False))
        f.write("\n\n")

rows = []
for gse, sub in df.groupby("GSE", sort=True):
    c = sub["Sample_review_status_v1"].value_counts().to_dict()
    rows.append({
        "GSE": gse,
        "N_GSM": sub["GSM"].nunique(),
        "N_HFD_CANDIDATE": c.get("HFD_CANDIDATE",0),
        "N_CONTROL_CANDIDATE": c.get("CONTROL_CANDIDATE",0),
        "N_REVIEW_WHOLE_LIVER": c.get("REVIEW_WHOLE_LIVER",0),
        "N_REVIEW_SPECIAL_DIET": c.get("REVIEW_SPECIAL_DIET",0),
        "N_REVIEW_TISSUE": c.get("REVIEW_TISSUE",0),
        "N_REVIEW_METADATA": c.get("REVIEW_METADATA",0),
        "N_RAW_UNAVAILABLE_SAMPLE": c.get("RAW_UNAVAILABLE_SAMPLE",0),
        "N_RAW_UNAVAILABLE_DATASET": c.get("RAW_UNAVAILABLE_DATASET",0),
    })

pd.DataFrame(rows).to_csv(SUMMARY, sep="\t", index=False)
print("Created:")
print(OUT)
print(SUMMARY)
