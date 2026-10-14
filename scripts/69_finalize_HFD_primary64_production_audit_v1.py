#!/usr/bin/env python3
from pathlib import Path
import hashlib
import pandas as pd
import numpy as np

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

MANIFEST = ROOT / "00_metadata/technical_harmonization/frozen_v1/HFD_primary_ENA_R1R2_download_manifest_frozen_v1.tsv"
PRIMARY_STATUS = ROOT / "00_metadata/processing_production/HFD_primary64_production_status_v1.tsv"
ADAPTER_RES = ROOT / "00_metadata/processing_production/exception_review/HFD_adapter_exception_downstream_status_v1.tsv"
STRAND_RES = ROOT / "00_metadata/processing_production/exception_review/HFD_strandedness_exception_resolution_status_v1.tsv"

OUTDIR = ROOT / "00_metadata/processing_production/final_audit_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_STATUS = OUTDIR / "HFD_primary64_final_resolved_status_v1.tsv"
OUT_GSE = OUTDIR / "HFD_primary64_final_matrix_audit_v1.tsv"
OUT_SUMMARY = OUTDIR / "HFD_primary64_final_production_audit_summary_v1.txt"

for q in [MANIFEST, PRIMARY_STATUS, ADAPTER_RES, STRAND_RES]:
    if not q.exists():
        raise RuntimeError(f"Missing required input: {q}")

m = pd.read_csv(MANIFEST, sep="\t", dtype=str).fillna("")
s = pd.read_csv(PRIMARY_STATUS, sep="\t", dtype=str).fillna("")
a = pd.read_csv(ADAPTER_RES, sep="\t", dtype=str).fillna("")
r = pd.read_csv(STRAND_RES, sep="\t", dtype=str).fillna("")

assert m["GSE"].nunique() == 64
assert m["GSM"].nunique() == 837
assert m["SRR"].nunique() == 951
assert len(m) == 1902

if set(m["FASTQ_component"]) != {"R1", "R2"}:
    raise RuntimeError(f"Unexpected FASTQ components: {set(m['FASTQ_component'])}")

per_srr_fastq = m.groupby("SRR")["FASTQ_component"].nunique()
if not (per_srr_fastq == 2).all():
    raise RuntimeError("At least one frozen SRR is not exactly paired R1/R2")

if s["GSE"].nunique() != 64 or len(s) != 64:
    raise RuntimeError("Primary production status must contain exactly 64 unique GSE")

adapter_set = set(a.loc[a["Resolution_status"].eq("PASS_AFTER_ADAPTER_TRIM"), "GSE"])
strand_set = set(r.loc[r["Resolution_status"].eq("PASS_AFTER_STRANDEDNESS_ADJUDICATION"), "GSE"])

if len(adapter_set) != 7:
    raise RuntimeError(f"Expected 7 adapter-resolved GSE, observed {len(adapter_set)}")
if len(strand_set) != 2:
    raise RuntimeError(f"Expected 2 strandedness-resolved GSE, observed {len(strand_set)}")
if adapter_set & strand_set:
    raise RuntimeError("Adapter and strandedness resolution sets overlap")

rows = []
for z in s.itertuples(index=False):
    gse = z.GSE
    original = z.Production_status
    if original == "PASS_OR_ALREADY_COMPLETE":
        route = "DIRECT_PRODUCTION_PASS"
        final = "PASS"
    elif original == "REVIEW_REQUIRED_ADAPTER":
        if gse not in adapter_set:
            raise RuntimeError(f"{gse}: adapter exception not resolved")
        route = "PASS_AFTER_ADAPTER_TRIM"
        final = "PASS"
    elif original == "REVIEW_REQUIRED_AMBIGUOUS_STRANDEDNESS":
        if gse not in strand_set:
            raise RuntimeError(f"{gse}: strandedness exception not resolved")
        route = "PASS_AFTER_STRANDEDNESS_ADJUDICATION"
        final = "PASS"
    else:
        raise RuntimeError(f"{gse}: unsupported unresolved original status {original}")
    rows.append({
        "GSE": gse,
        "Original_production_status": original,
        "Final_resolution_route": route,
        "Final_status": final,
    })

final_status = pd.DataFrame(rows).sort_values("GSE")
if len(final_status) != 64 or final_status["GSE"].nunique() != 64:
    raise RuntimeError("Final resolved status is not 64 unique GSE")
if not (final_status["Final_status"] == "PASS").all():
    raise RuntimeError("Not all final statuses are PASS")
final_status.to_csv(OUT_STATUS, sep="\t", index=False)

audit_rows = []
reference_gene_hash = None
reference_gene_n = None

for gse in sorted(m["GSE"].unique()):
    expected_gsm = sorted(m.loc[m["GSE"].eq(gse), "GSM"].unique())
    expected_srr = sorted(m.loc[m["GSE"].eq(gse), "SRR"].unique())

    matrix = ROOT / "03_counts/gsm" / gse / f"{gse}_GSM_counts_v1.tsv"
    sample_map = ROOT / "03_counts/gsm" / gse / f"{gse}_sample_map_v1.tsv"
    done = ROOT / "00_metadata/processing_production" / gse / f"{gse}_PRODUCTION_DONE_v1.txt"

    if not matrix.exists():
        raise RuntimeError(f"{gse}: missing final GSM matrix")
    if not sample_map.exists():
        raise RuntimeError(f"{gse}: missing final sample map")
    if not done.exists():
        raise RuntimeError(f"{gse}: missing PRODUCTION_DONE marker")

    x = pd.read_csv(matrix, sep="\t")
    if "Geneid" not in x.columns:
        raise RuntimeError(f"{gse}: matrix lacks Geneid")
    if x["Geneid"].duplicated().any():
        raise RuntimeError(f"{gse}: duplicate Geneid")

    observed_gsm = sorted([c for c in x.columns if c != "Geneid"])
    missing_gsm = sorted(set(expected_gsm) - set(observed_gsm))
    extra_gsm = sorted(set(observed_gsm) - set(expected_gsm))
    if missing_gsm or extra_gsm:
        raise RuntimeError(f"{gse}: GSM mismatch missing={missing_gsm} extra={extra_gsm}")

    vals = x[observed_gsm].apply(pd.to_numeric, errors="raise")
    arr = vals.to_numpy(dtype=float)
    if not np.isfinite(arr).all():
        raise RuntimeError(f"{gse}: non-finite counts")
    if (arr < 0).any():
        raise RuntimeError(f"{gse}: negative counts")
    if not np.equal(arr, np.floor(arr)).all():
        raise RuntimeError(f"{gse}: non-integer counts")

    sm = pd.read_csv(sample_map, sep="\t", dtype=str).fillna("")
    observed_srr = sorted(sm["SRR"].unique())
    observed_map_gsm = sorted(sm["GSM"].unique())

    if observed_srr != expected_srr:
        raise RuntimeError(f"{gse}: sample-map SRR set differs from frozen manifest")
    if observed_map_gsm != expected_gsm:
        raise RuntimeError(f"{gse}: sample-map GSM set differs from frozen manifest")

    exp_map = (
        m.loc[m["GSE"].eq(gse), ["GSM","SRR","Contrast_id","Role"]]
        .drop_duplicates()
        .sort_values(["GSM","SRR","Contrast_id","Role"])
        .reset_index(drop=True)
    )
    obs_map = (
        sm[["GSM","SRR","Contrast_id","Role"]]
        .drop_duplicates()
        .sort_values(["GSM","SRR","Contrast_id","Role"])
        .reset_index(drop=True)
    )
    if not exp_map.equals(obs_map):
        raise RuntimeError(f"{gse}: final sample map differs from frozen biological mapping")

    gene_bytes = ("\n".join(x["Geneid"].astype(str)) + "\n").encode()
    gene_hash = hashlib.sha256(gene_bytes).hexdigest()
    gene_n = len(x)

    if reference_gene_hash is None:
        reference_gene_hash = gene_hash
        reference_gene_n = gene_n
    else:
        if gene_hash != reference_gene_hash:
            raise RuntimeError(f"{gse}: Geneid/order differs across final matrices")
        if gene_n != reference_gene_n:
            raise RuntimeError(f"{gse}: gene-row count differs across final matrices")

    route = final_status.loc[
        final_status["GSE"].eq(gse), "Final_resolution_route"
    ].iloc[0]

    audit_rows.append({
        "GSE": gse,
        "Final_resolution_route": route,
        "Expected_GSM": len(expected_gsm),
        "Observed_GSM": len(observed_gsm),
        "Expected_SRR": len(expected_srr),
        "Observed_SRR": len(observed_srr),
        "Gene_rows": gene_n,
        "Geneid_SHA256": gene_hash,
        "DONE_marker": "PASS",
        "Matrix_QC": "PASS",
        "Frozen_mapping_QC": "PASS",
    })

audit = pd.DataFrame(audit_rows).sort_values("GSE")
audit.to_csv(OUT_GSE, sep="\t", index=False)

done_files = list((ROOT / "00_metadata/processing_production").glob("*/GSE*_PRODUCTION_DONE_v1.txt"))
matrix_files = list((ROOT / "03_counts/gsm").glob("GSE*/GSE*_GSM_counts_v1.tsv"))

if len(done_files) != 64:
    raise RuntimeError(f"Expected 64 DONE markers, observed {len(done_files)}")
if len(matrix_files) != 64:
    raise RuntimeError(f"Expected 64 final GSM matrices, observed {len(matrix_files)}")
if audit["Observed_GSM"].sum() != 837:
    raise RuntimeError(f"Expected 837 GSM across matrices, observed {audit['Observed_GSM'].sum()}")
if audit["Observed_SRR"].sum() != 951:
    raise RuntimeError(f"Expected 951 SRR across sample maps, observed {audit['Observed_SRR'].sum()}")

route_counts = final_status["Final_resolution_route"].value_counts().to_dict()

summary = [
    "===== HFD PRIMARY64 FINAL PRODUCTION AUDIT v1 =====",
    "",
    "FINAL_STATUS=PASS",
    f"Independent_GSE={len(audit)}",
    f"Frozen_GSM={m['GSM'].nunique()}",
    f"Frozen_SRR={m['SRR'].nunique()}",
    f"Frozen_FASTQ={len(m)}",
    f"Final_DONE_markers={len(done_files)}",
    f"Final_GSM_matrices={len(matrix_files)}",
    f"Observed_GSM_across_final_matrices={int(audit['Observed_GSM'].sum())}",
    f"Observed_SRR_across_final_sample_maps={int(audit['Observed_SRR'].sum())}",
    f"Common_gene_rows={reference_gene_n}",
    f"Common_Geneid_SHA256={reference_gene_hash}",
    "",
    f"DIRECT_PRODUCTION_PASS={route_counts.get('DIRECT_PRODUCTION_PASS', 0)}",
    f"PASS_AFTER_ADAPTER_TRIM={route_counts.get('PASS_AFTER_ADAPTER_TRIM', 0)}",
    f"PASS_AFTER_STRANDEDNESS_ADJUDICATION={route_counts.get('PASS_AFTER_STRANDEDNESS_ADJUDICATION', 0)}",
    "",
    "MATRIX_QC=64/64 PASS",
    "FROZEN_GSM_MEMBERSHIP_QC=64/64 PASS",
    "FROZEN_SRR_MAPPING_QC=64/64 PASS",
    "GENEID_ORDER_CONSISTENCY=64/64 PASS",
    "PRIMARY64_MASTER_STATUS_V1_MODIFIED=NO",
    "",
    f"Final_resolved_status={OUT_STATUS.relative_to(ROOT)}",
    f"Final_matrix_audit={OUT_GSE.relative_to(ROOT)}",
]

OUT_SUMMARY.write_text("\n".join(summary) + "\n")
print("\n".join(summary))

