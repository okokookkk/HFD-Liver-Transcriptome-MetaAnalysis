#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()
F = ROOT / "00_metadata/technical_harmonization/frozen_v1/HFD_primary_ENA_R1R2_download_manifest_frozen_v1.tsv"
OUT = ROOT / "00_metadata/processing_production"
OUT.mkdir(parents=True, exist_ok=True)

x = pd.read_csv(F, sep="\t", dtype=str).fillna("")
x["FASTQ_bytes"] = pd.to_numeric(x["FASTQ_bytes"], errors="raise")

q = (
    x.groupby("GSE", as_index=False)
     .agg(
         GSM=("GSM", "nunique"),
         SRR=("SRR", "nunique"),
         Contrasts=("Contrast_id", "nunique"),
         FASTQ_files=("FASTQ_filename", "nunique"),
         FASTQ_bytes=("FASTQ_bytes", "sum"),
     )
)

q["FASTQ_GiB"] = q["FASTQ_bytes"] / 1024**3
q["Pilot_status"] = q["GSE"].map({
    "GSE268200": "PILOT1_VALIDATED",
    "GSE151750": "PILOT2_VALIDATED",
}).fillna("")

q = q.sort_values(["FASTQ_bytes", "GSE"]).reset_index(drop=True)
q.insert(0, "Queue_order", range(1, len(q) + 1))

assert q["GSE"].nunique() == 64
assert q["GSM"].sum() == 837
assert q["SRR"].sum() == 951
assert q["FASTQ_files"].sum() == 1902

path = OUT / "HFD_primary64_production_queue_v1.tsv"
q.to_csv(path, sep="\t", index=False)

print("===== PRODUCTION QUEUE v1 =====")
print("GSE:", q["GSE"].nunique())
print("GSM:", q["GSM"].sum())
print("SRR:", q["SRR"].sum())
print("FASTQ files:", q["FASTQ_files"].sum())
print("FASTQ TiB:", q["FASTQ_bytes"].sum() / 1024**4)
print("\nFirst 15:")
print(q.head(15).to_string(index=False))
print("\nSaved:", path)
