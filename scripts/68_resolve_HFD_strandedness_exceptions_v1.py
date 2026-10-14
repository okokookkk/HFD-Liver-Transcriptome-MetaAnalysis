#!/usr/bin/env python3
from pathlib import Path
from datetime import datetime
import argparse
import shutil
import subprocess
import pandas as pd
import numpy as np

ROOT_DEFAULT = __import__("os").environ.get("HFD_PROJECT_ROOT", ".")
TARGETS = ["GSE137022", "GSE278920"]

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=ROOT_DEFAULT)
ap.add_argument("--fc-threads", type=int, default=4)
args = ap.parse_args()

ROOT = Path(args.root)

MANIFEST = ROOT / "00_metadata/technical_harmonization/frozen_v1/HFD_primary_ENA_R1R2_download_manifest_frozen_v1.tsv"
PRIMARY_STATUS = ROOT / "00_metadata/processing_production/HFD_primary64_production_status_v1.tsv"
EXC_ROOT = ROOT / "00_metadata/processing_production/exception_review"
GTF = ROOT / "00_metadata/reference/GENCODE_M39_GRCm39/gencode.vM39.primary_assembly.annotation.gtf"
FC = Path(__import__('os').environ.get('HFD_FEATURECOUNTS', 'featureCounts'))

for q in [MANIFEST, PRIMARY_STATUS, GTF, FC]:
    if not q.exists():
        raise RuntimeError(f"Required input/tool not found: {q}")

manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str).fillna("")
status = pd.read_csv(PRIMARY_STATUS, sep="\t", dtype=str).fillna("")

assert manifest["GSE"].nunique() == 64
assert status["GSE"].nunique() == 64

def now():
    return datetime.now().astimezone().isoformat()

def parse_fc_counts(path, srr):
    x = pd.read_csv(path, sep="\t", comment="#")
    meta_cols = {"Geneid","Chr","Start","End","Strand","Length"}
    count_cols = [c for c in x.columns if c not in meta_cols]
    if len(count_cols) != 1:
        raise RuntimeError(f"{srr}: expected exactly one count column; got {count_cols}")
    out = x[["Geneid", count_cols[0]]].copy()
    out = out.rename(columns={count_cols[0]: srr})
    out[srr] = pd.to_numeric(out[srr], errors="raise").astype(np.int64)
    if (out[srr] < 0).any():
        raise RuntimeError(f"{srr}: negative counts")
    return out

def parse_fc_summary(path):
    x = pd.read_csv(path, sep="\t", index_col=0)
    if x.shape[1] != 1:
        raise RuntimeError(f"{path}: expected one summary column")
    return pd.to_numeric(x.iloc[:,0], errors="raise")

audit_rows = []
batch_rows = []

for gse in TARGETS:
    print("\n" + "=" * 100)
    print(f"{gse}: STRANDEDNESS EXCEPTION RESOLUTION START")

    z = status[status["GSE"].eq(gse)]
    if len(z) != 1 or z.iloc[0]["Production_status"] != "REVIEW_REQUIRED_AMBIGUOUS_STRANDEDNESS":
        raise RuntimeError(f"{gse}: primary status is not REVIEW_REQUIRED_AMBIGUOUS_STRANDEDNESS")

    mx = manifest[manifest["GSE"].eq(gse)].copy()
    if mx.empty:
        raise RuntimeError(f"{gse}: absent from frozen manifest")

    runmap = (
        mx[["GSE","GSM","SRR","Contrast_id","Role"]]
        .drop_duplicates()
        .sort_values(["GSM","SRR"])
    )

    if runmap.groupby("SRR")["GSM"].nunique().max() != 1:
        raise RuntimeError(f"{gse}: SRR maps to multiple GSM")
    if runmap.groupby("GSM")["Role"].nunique().max() != 1:
        raise RuntimeError(f"{gse}: GSM maps to multiple Role")
    if runmap.groupby("GSM")["Contrast_id"].nunique().max() != 1:
        raise RuntimeError(f"{gse}: GSM maps to multiple Contrast_id")

    expected_srr = sorted(runmap["SRR"].unique())
    expected_gsm = sorted(runmap["GSM"].unique())

    star_dir = ROOT / "03_counts/star" / gse
    starqc_dir = ROOT / "02_qc/star" / gse
    fc_dir = ROOT / "03_counts/featurecounts" / gse / "per_srr_strand_adjudicated_v1"
    gsm_dir = ROOT / "03_counts/gsm" / gse
    log_dir = ROOT / "logs/production_exception_review" / gse
    prod_meta = ROOT / "00_metadata/processing_production" / gse
    exc_meta = EXC_ROOT / gse

    for d in [fc_dir, gsm_dir, log_dir, prod_meta, exc_meta]:
        d.mkdir(parents=True, exist_ok=True)

    empirical_file = starqc_dir / f"{gse}_STAR_strandedness_v1.tsv"
    if not empirical_file.exists():
        raise RuntimeError(f"{gse}: missing empirical strandedness table")

    empirical = pd.read_csv(empirical_file, sep="\t")
    if set(empirical["SRR"]) != set(expected_srr):
        raise RuntimeError(f"{gse}: empirical SRR set does not match frozen manifest")

    # Study-level manual adjudication:
    # Both publications explicitly report Illumina TruSeq Stranded library preparation.
    # For standard paired-end TruSeq Stranded RNA-seq, featureCounts orientation is -s 2.
    adjud = empirical.copy()
    adjud["Empirical_featureCounts_s"] = adjud["featureCounts_s"]
    adjud["Empirical_Strand_call"] = adjud["Strand_call"]
    adjud["Adjudicated_featureCounts_s"] = 2
    adjud["Adjudicated_Strand_call"] = "REVERSE_STRANDED"
    adjud["Adjudication_level"] = "GSE_PROTOCOL"
    adjud["Adjudication_reason"] = (
        "Publication explicitly reports Illumina TruSeq Stranded library preparation; "
        "use one study-consistent reverse-stranded counting mode for all frozen SRRs."
    )

    adjud_file = exc_meta / f"{gse}_strandedness_manual_adjudication_v1.tsv"
    adjud.to_csv(adjud_file, sep="\t", index=False)

    for r in adjud.itertuples(index=False):
        audit_rows.append({
            "GSE": gse,
            "GSM": r.GSM,
            "SRR": r.SRR,
            "Forward_fraction": r.Forward_fraction,
            "Reverse_fraction": r.Reverse_fraction,
            "Empirical_Strand_call": r.Empirical_Strand_call,
            "Empirical_featureCounts_s": r.Empirical_featureCounts_s,
            "Adjudicated_Strand_call": r.Adjudicated_Strand_call,
            "Adjudicated_featureCounts_s": r.Adjudicated_featureCounts_s,
            "Adjudication_level": r.Adjudication_level,
            "Adjudication_reason": r.Adjudication_reason,
        })

    # Recompute featureCounts for ALL SRRs in the GSE using s=2.
    srr_tables = []
    fc_qc_rows = []

    for srr in expected_srr:
        bam = star_dir / srr / "Aligned.sortedByCoord.out.bam"
        if not bam.exists():
            raise RuntimeError(f"{gse}/{srr}: missing STAR BAM")

        out = fc_dir / f"{srr}_featureCounts_strand_adjudicated_v1.txt"
        summary = Path(str(out) + ".summary")
        log = log_dir / f"{srr}_featureCounts_strand_adjudicated_v1.log"

        # Always recompute in this dedicated adjudicated directory.
        cmd = [
            str(FC),
            "-T", str(args.fc_threads),
            "-p", "--countReadPairs", "-B", "-C",
            "-s", "2",
            "-t", "exon",
            "-g", "gene_id",
            "-a", str(GTF),
            "-o", str(out),
            str(bam),
        ]

        print(f"{srr}: featureCounts -s 2")

        with log.open("w") as fh:
            rc = subprocess.run(
                cmd,
                stdout=fh,
                stderr=subprocess.STDOUT
            ).returncode

        if rc != 0 or not (out.exists() and summary.exists()):
            raise RuntimeError(f"{gse}/{srr}: featureCounts failed; see {log}")

        c = parse_fc_counts(out, srr)
        srr_tables.append(c)

        sm = parse_fc_summary(summary)
        assigned = int(sm.get("Assigned", 0))
        total = int(sm.sum())

        fc_qc_rows.append({
            "GSE": gse,
            "SRR": srr,
            "featureCounts_s": 2,
            "Assigned_fragments": assigned,
            "FC_summary_total_alignments": total,
            "FC_Unassigned_MultiMapping_alignments": int(sm.get("Unassigned_MultiMapping", 0)),
            "FC_Unassigned_NoFeatures": int(sm.get("Unassigned_NoFeatures", 0)),
            "FC_Unassigned_Ambiguity": int(sm.get("Unassigned_Ambiguity", 0)),
        })

    fcqc = pd.DataFrame(fc_qc_rows).merge(
        runmap[["GSM","SRR","Contrast_id","Role"]],
        on="SRR",
        how="left"
    ).sort_values(["GSM","SRR"])

    fcqc_file = starqc_dir / f"{gse}_featureCounts_QC_strand_adjudicated_v1.tsv"
    fcqc.to_csv(fcqc_file, sep="\t", index=False)

    # Merge SRR counts by exact Geneid.
    srr_counts = srr_tables[0]
    for q in srr_tables[1:]:
        srr_counts = srr_counts.merge(q, on="Geneid", how="inner", validate="one_to_one")

    expected_cols = ["Geneid"] + expected_srr
    if set(srr_counts.columns) != set(expected_cols):
        raise RuntimeError(f"{gse}: SRR matrix columns do not match frozen SRR set")
    srr_counts = srr_counts[expected_cols]

    # Aggregate technical SRRs to biological GSM by exact summation.
    gsm_counts = pd.DataFrame({"Geneid": srr_counts["Geneid"]})
    agg_rows = []

    for gsm, z2 in runmap.groupby("GSM", sort=True):
        srrs = sorted(z2["SRR"].tolist())
        role = z2["Role"].unique()
        contrast = z2["Contrast_id"].unique()

        if len(role) != 1 or len(contrast) != 1:
            raise RuntimeError(f"{gse}/{gsm}: frozen GSM mapping inconsistency")

        gsm_counts[gsm] = srr_counts[srrs].sum(axis=1).astype(np.int64)

        before = int(srr_counts[srrs].to_numpy().sum())
        after = int(gsm_counts[gsm].sum())

        agg_rows.append({
            "GSM": gsm,
            "Role": role[0],
            "Contrast_id": contrast[0],
            "SRR_count": len(srrs),
            "SRRs": ",".join(srrs),
            "SRR_total_assigned": before,
            "GSM_total_after_sum": after,
            "Difference": after - before,
            "Aggregation_QC": "PASS" if after == before else "FAIL",
        })

    agg = pd.DataFrame(agg_rows).sort_values("GSM")
    agg_file = starqc_dir / f"{gse}_SRR_to_GSM_aggregation_QC_strand_adjudicated_v1.tsv"
    agg.to_csv(agg_file, sep="\t", index=False)

    if not (agg["Aggregation_QC"] == "PASS").all():
        raise RuntimeError(f"{gse}: SRR-to-GSM aggregation failure")

    # Canonical final GSM matrix and sample map.
    final = gsm_dir / f"{gse}_GSM_counts_v1.tsv"
    sample_map = gsm_dir / f"{gse}_sample_map_v1.tsv"

    gsm_counts = gsm_counts[["Geneid"] + expected_gsm]
    gsm_counts.to_csv(final, sep="\t", index=False)
    runmap.to_csv(sample_map, sep="\t", index=False)

    if gsm_counts["Geneid"].duplicated().any():
        raise RuntimeError(f"{gse}: duplicate Geneid in final matrix")
    if (gsm_counts.drop(columns="Geneid") < 0).any().any():
        raise RuntimeError(f"{gse}: negative final counts")

    exc_status = exc_meta / f"{gse}_STRANDEDNESS_EXCEPTION_RESOLVED_v1.txt"
    exc_status.write_text(
        "\n".join([
            f"GSE={gse}",
            "STATUS=PASS_AFTER_STRANDEDNESS_ADJUDICATION",
            f"TIMESTAMP={now()}",
            "ADJUDICATED_FEATURECOUNTS_S=2",
            "ADJUDICATION_LEVEL=GSE_PROTOCOL",
            f"SRR={len(expected_srr)}",
            f"GSM={len(expected_gsm)}",
            f"MANUAL_ADJUDICATION_TABLE={adjud_file.relative_to(ROOT)}",
            f"FEATURECOUNTS_QC={fcqc_file.relative_to(ROOT)}",
            f"AGGREGATION_QC={agg_file.relative_to(ROOT)}",
            f"FINAL={final.relative_to(ROOT)}",
            "RAW_FASTQ_MODIFIED=NO",
            "STAR_ALIGNMENT_REUSED=YES",
            "FROZEN_TECHNICAL_ELIGIBILITY_CHANGED=NO",
            "FROZEN_GSM_SRR_CONTRAST_MEMBERSHIP_CHANGED=NO",
            "PRIMARY64_MASTER_STATUS_V1_MODIFIED=NO",
        ]) + "\n"
    )

    done = prod_meta / f"{gse}_PRODUCTION_DONE_v1.txt"
    done.write_text(
        "\n".join([
            f"GSE={gse}",
            "STATUS=PASS",
            f"TIMESTAMP={now()}",
            "RESOLVED_FROM=REVIEW_REQUIRED_AMBIGUOUS_STRANDEDNESS",
            "RESOLUTION=STUDY_LEVEL_TRUSEQ_STRANDED_PROTOCOL_ADJUDICATION_S2",
            "DOWNSTREAM_DRIVER=68_resolve_HFD_strandedness_exceptions_v1.py",
            f"FINAL={final.relative_to(ROOT)}",
            f"SAMPLE_MAP={sample_map.relative_to(ROOT)}",
            f"STRANDEDNESS_ADJUDICATION={adjud_file.relative_to(ROOT)}",
            f"FEATURECOUNTS_QC={fcqc_file.relative_to(ROOT)}",
            f"AGGREGATION_QC={agg_file.relative_to(ROOT)}",
            "RAW_FASTQ_MODIFIED=NO",
            "STAR_ALIGNMENT_REUSED=YES",
            "FROZEN_TECHNICAL_ELIGIBILITY_CHANGED=NO",
        ]) + "\n"
    )

    batch_rows.append({
        "GSE": gse,
        "Exit_code": 0,
        "Resolution_status": "PASS_AFTER_STRANDEDNESS_ADJUDICATION",
        "Adjudicated_featureCounts_s": 2,
        "SRR": len(expected_srr),
        "GSM": len(expected_gsm),
    })

    print(f"===== {gse} PASS AFTER STRANDEDNESS ADJUDICATION =====")
    print("GSM matrix:", final)

audit = pd.DataFrame(audit_rows)
audit_out = EXC_ROOT / "HFD_strandedness_manual_adjudication_v1.tsv"
audit.to_csv(audit_out, sep="\t", index=False)

batch = pd.DataFrame(batch_rows)
batch_out = EXC_ROOT / "HFD_strandedness_exception_resolution_status_v1.tsv"
batch.to_csv(batch_out, sep="\t", index=False)

print("\n" + "=" * 100)
print("STRANDEDNESS EXCEPTION RESOLUTION COMPLETE")
print(batch.to_string(index=False))
print("\nSaved:")
print(" ", audit_out.relative_to(ROOT))
print(" ", batch_out.relative_to(ROOT))
print("Primary64 master status v1 was NOT modified.")

