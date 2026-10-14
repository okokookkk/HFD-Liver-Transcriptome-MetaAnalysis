#!/usr/bin/env python3
"""
04_repair_ssRNAseq_and_rebuild_acquisition_v1_1.py

Administrative/raw-phase parser correction:
- v1 recognized only Library_strategy == "RNA-Seq".
- GSE296500 contains 45 whole-liver PolyA RNA-seq samples labeled "ssRNA-seq".
- This script rescues their SRX->SRR mappings, creates a new run table v1.1,
  and rebuilds dissertation-grade acquisition/review outputs without modifying v1.

v1 remains immutable for provenance.
"""

from pathlib import Path
from datetime import datetime
import csv
import io
import hashlib
import re
import subprocess
import time
import pandas as pd

ROOT = Path(".")
MAP = ROOT / "00_metadata" / "run_mapping"
FROZEN = ROOT / "00_metadata" / "frozen_search"
ACQ = ROOT / "00_metadata" / "data_acquisition"
SCREEN = ROOT / "00_metadata" / "sample_screening"
LOGDIR = ROOT / "logs"

for p in [ACQ, SCREEN, LOGDIR]:
    p.mkdir(parents=True, exist_ok=True)

GEO = MAP / "primary107_geo_samples_v1.tsv"
RUNS_V1 = MAP / "primary107_sra_runs_v1.tsv"
EXC = MAP / "primary107_raw_availability_exceptions_v1.tsv"
DEC = MAP / "primary107_raw_availability_decisions_v1.tsv"
PRIMARY = FROZEN / "HFD_primary_analysis_GSE_frozen_v1.txt"

RUNS_V11 = MAP / "primary107_sra_runs_v1_1.tsv"
RESCUED = MAP / "primary107_ssRNAseq_rescued_runs_v1_1.tsv"
UNRESOLVED = MAP / "primary107_ssRNAseq_unresolved_v1_1.tsv"
MAP_SUMMARY = MAP / "primary107_mapping_summary_v1_1.txt"
CORRECTION = MAP / "DATA_ACQUISITION_CORRECTION_v1.1.txt"

MASTER = ACQ / "Data_Acquisition_Master_v1_1.tsv"
RUN_MASTER = ACQ / "Data_Acquisition_Run_Master_v1_1.tsv"
GSE_SUMMARY = ACQ / "Data_Acquisition_GSE_Summary_v1_1.tsv"
REVIEW_QUEUE = ACQ / "Data_Acquisition_Review_Queue_v1_1.tsv"
MANIFEST = ACQ / "Data_Acquisition_Manifest_v1_1.txt"

REVIEW_PACK = SCREEN / "primary107_sample_review_pack_v1_1.txt"
REVIEW_SUMMARY = SCREEN / "primary107_sample_review_summary_v1_1.tsv"

NOW = datetime.now().astimezone().isoformat()
BULK_RNA_STRATEGIES = {"RNA-SEQ", "SSRNA-SEQ"}


def norm(x):
    return "" if pd.isna(x) else str(x).strip()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def uniq_join(series):
    return ";".join(sorted(set(norm(x) for x in series if norm(x))))


def is_bulk_rna(x):
    return norm(x).upper() in BULK_RNA_STRATEGIES


def fetch_runinfo(srx):
    cmd = f'esearch -db sra -query "{srx}" | efetch -format runinfo'
    p = subprocess.run(
        ["bash", "-lc", cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=180,
    )
    if p.returncode != 0:
        return [], p.stderr.strip()

    text = p.stdout.strip()
    if not text:
        return [], ""

    try:
        rows = list(csv.DictReader(io.StringIO(text)))
    except Exception as e:
        return [], f"CSV parse error: {e}"

    exact = [
        r for r in rows
        if norm(r.get("Experiment", "")) == srx
    ]
    return exact, ""


def contains_any(text, patterns):
    t = norm(text).lower()
    return any(re.search(p, t, flags=re.I) for p in patterns)


geo = pd.read_csv(GEO, sep="\t", dtype=str).fillna("")
runs_v1 = pd.read_csv(RUNS_V1, sep="\t", dtype=str).fillna("")
primary = [x.strip() for x in PRIMARY.read_text().splitlines() if x.strip()]

assert len(primary) == 107
assert geo["GSE"].nunique() == 107

# ------------------------------------------------------------------
# Identify parser-missed bulk RNA-seq samples
# ------------------------------------------------------------------
bulk = geo[geo["Library_strategy"].map(is_bulk_rna)].copy()
nonstandard = bulk[bulk["Library_strategy"].str.upper().ne("RNA-SEQ")].copy()

mapped_gsm_v1 = set(runs_v1["GSM"])
need_rescue = nonstandard[
    (~nonstandard["GSM"].isin(mapped_gsm_v1)) &
    (nonstandard["SRX"].str.strip().ne(""))
].copy()

print("Bulk RNA-seq GSM (RNA-Seq + ssRNA-seq):", len(bulk))
print("Nonstandard bulk RNA-seq GSM:", len(nonstandard))
print("Need SRX->SRR rescue:", len(need_rescue))
print("Nonstandard GSE:", ",".join(sorted(nonstandard["GSE"].unique())))

# Current known correction should be isolated to GSE296500 / 45 GSM.
if len(nonstandard) != 45 or set(nonstandard["GSE"]) != {"GSE296500"}:
    raise SystemExit(
        "Unexpected nonstandard bulk RNA-seq universe. "
        "Review Library_strategy values before proceeding."
    )

# ------------------------------------------------------------------
# Rescue SRX->SRR
# ------------------------------------------------------------------
rescued_rows = []
unresolved_rows = []

existing_cols = list(runs_v1.columns)

for idx, (_, r) in enumerate(need_rescue.iterrows(), start=1):
    gse, gsm, srx = r["GSE"], r["GSM"], r["SRX"]
    print(f"[{idx:02d}/{len(need_rescue)}] {gse} {gsm} {srx}")

    rows, err = fetch_runinfo(srx)

    if not rows:
        unresolved_rows.append({
            "GSE": gse,
            "GSM": gsm,
            "SRX": srx,
            "Library_strategy": r["Library_strategy"],
            "Status": "UNRESOLVED",
            "Error": err,
        })
        print("   UNRESOLVED", err)
        time.sleep(0.35)
        continue

    for rr in rows:
        out = {c: "" for c in existing_cols}
        mapping = {
            "GSE": gse,
            "GSM": gsm,
            "SRX": norm(rr.get("Experiment", "")),
            "SRR": norm(rr.get("Run", "")),
            "BioSample": norm(rr.get("BioSample", "")),
            "BioProject": norm(rr.get("BioProject", "")),
            "LibraryLayout": norm(rr.get("LibraryLayout", "")),
            "LibraryStrategy": norm(rr.get("LibraryStrategy", "")),
            "LibrarySource": norm(rr.get("LibrarySource", "")),
            "Platform": norm(rr.get("Platform", "")),
            "Model": norm(rr.get("Model", "")),
            "spots": norm(rr.get("spots", "")),
            "bases": norm(rr.get("bases", "")),
            "avgLength": norm(rr.get("avgLength", "")),
            "size_MB": norm(rr.get("size_MB", "")),
            "ReleaseDate": norm(rr.get("ReleaseDate", "")),
            "Sample_title": r.get("Sample_title", ""),
            "Source_name": r.get("Source_name", ""),
            "Tissue": r.get("Tissue", ""),
            "Diet": r.get("Diet", ""),
            "Genotype": r.get("Genotype", ""),
            "Sex": r.get("Sex", ""),
            "Age": r.get("Age", ""),
            "Strain": r.get("Strain", ""),
            "Treatment": r.get("Treatment", ""),
        }
        for k, v in mapping.items():
            if k in out:
                out[k] = v
        rescued_rows.append(out)

    print("   RESCUED runs:", len(rows))
    time.sleep(0.35)

rescued_df = pd.DataFrame(rescued_rows, columns=existing_cols)
unresolved_df = pd.DataFrame(unresolved_rows)

rescued_df.to_csv(RESCUED, sep="\t", index=False)
unresolved_df.to_csv(UNRESOLVED, sep="\t", index=False)

if len(unresolved_df) > 0:
    raise SystemExit(
        f"{len(unresolved_df)} ssRNA-seq GSM remain unresolved. "
        f"See {UNRESOLVED}. Do not rebuild/freeze v1.1 yet."
    )

# ------------------------------------------------------------------
# Append + validate run table
# ------------------------------------------------------------------
runs = pd.concat([runs_v1, rescued_df], ignore_index=True).fillna("")

if runs["SRR"].duplicated().any():
    dup = runs[runs["SRR"].duplicated(keep=False)].sort_values("SRR")
    dup.to_csv(MAP / "primary107_srr_duplicate_error_v1_1.tsv", sep="\t", index=False)
    raise SystemExit("Duplicate SRR detected after rescue.")

runs.to_csv(RUNS_V11, sep="\t", index=False)

# ------------------------------------------------------------------
# Rebuild acquisition master v1.1
# ------------------------------------------------------------------
for c in [
    "SRX", "SRR", "BioSample", "BioProject", "LibraryLayout",
    "LibraryStrategy", "LibrarySource", "Platform", "Model", "ReleaseDate"
]:
    if c not in runs.columns:
        runs[c] = ""

agg = (
    runs.groupby(["GSE", "GSM"], dropna=False)
    .agg(
        N_SRX=("SRX", lambda s: len(set(x for x in s if norm(x)))),
        SRX_list=("SRX", uniq_join),
        N_SRR=("SRR", lambda s: len(set(x for x in s if norm(x)))),
        SRR_list=("SRR", uniq_join),
        BioSample_SRA=("BioSample", uniq_join),
        BioProject=("BioProject", uniq_join),
        LibraryLayout=("LibraryLayout", uniq_join),
        SRA_LibraryStrategy=("LibraryStrategy", uniq_join),
        SRA_LibrarySource=("LibrarySource", uniq_join),
        Platform=("Platform", uniq_join),
        Model=("Model", uniq_join),
        ReleaseDate=("ReleaseDate", uniq_join),
    )
    .reset_index()
)

m = geo.merge(agg, on=["GSE", "GSM"], how="left").fillna("")

for c in [
    "Raw_availability", "Raw_primary_status", "Raw_availability_reason",
    "Processed_data_status", "GSE_raw_harmonized_status"
]:
    m[c] = ""

if EXC.exists() and EXC.stat().st_size > 0:
    ex = pd.read_csv(EXC, sep="\t", dtype=str).fillna("")
    exmap = ex.set_index(["GSE", "GSM"]).to_dict("index")
    for i, r in m.iterrows():
        key = (r["GSE"], r["GSM"])
        if key in exmap:
            e = exmap[key]
            m.at[i, "Raw_availability"] = norm(e.get("Raw_availability", ""))
            m.at[i, "Raw_primary_status"] = norm(e.get("Raw_primary_status", ""))
            m.at[i, "Raw_availability_reason"] = norm(e.get("Reason", ""))
            m.at[i, "Processed_data_status"] = norm(e.get("Processed_data_status", ""))

for i, r in m.iterrows():
    if m.at[i, "Raw_availability"]:
        continue

    bulk_rna = is_bulk_rna(r.get("Library_strategy", ""))
    n_srr = int(float(norm(r.get("N_SRR", "")) or 0))

    if bulk_rna and n_srr > 0:
        m.at[i, "Raw_availability"] = "PUBLIC_SRA"
        m.at[i, "Raw_primary_status"] = "RAW_MAPPED"
    elif bulk_rna:
        m.at[i, "Raw_availability"] = "NO_DIRECT_PUBLIC_SRA"
        m.at[i, "Raw_primary_status"] = "REVIEW_REQUIRED"
    else:
        m.at[i, "Raw_availability"] = "NOT_APPLICABLE_NON_BULK_RNASEQ"
        m.at[i, "Raw_primary_status"] = "NOT_PRIMARY_RNASEQ"

if DEC.exists() and DEC.stat().st_size > 0:
    d = pd.read_csv(DEC, sep="\t", dtype=str).fillna("")
    dmap = dict(zip(d["GSE"], d["Raw_harmonized_status"]))
    m["GSE_raw_harmonized_status"] = m["GSE"].map(dmap).fillna("")

m.loc[
    m["GSE_raw_harmonized_status"].eq(""),
    "GSE_raw_harmonized_status"
] = "RAW_ANALYZABLE_PENDING_SAMPLE_REVIEW"

combo_cols = [
    "Sample_title", "Source_name", "Tissue", "Diet", "Genotype",
    "Treatment", "All_characteristics", "All_relations"
]
m["_combined"] = m.apply(
    lambda r: " || ".join(norm(r.get(c, "")) for c in combo_cols),
    axis=1
)

m["Flag_mouse"] = (
    m["Organism"].str.lower().eq("mus musculus")
    .map({True: "YES", False: "NO"})
)
m["Flag_bulk_RNAseq"] = (
    m["Library_strategy"].map(is_bulk_rna)
    .map({True: "YES", False: "NO"})
)
# keep legacy-compatible column too
m["Flag_RNAseq"] = m["Flag_bulk_RNAseq"]

m["Flag_liver_text"] = m["_combined"].apply(
    lambda x: "YES" if contains_any(x, [r"\bliver\b", r"hepatic"]) else "NO"
)
m["Flag_non_whole_liver_signal"] = m["_combined"].apply(
    lambda x: "YES" if contains_any(
        x,
        [
            r"hepatocyte", r"kupffer", r"stellate", r"endothelial",
            r"cholangiocyte", r"single[- ]cell", r"single[- ]nucleus",
            r"tumou?r", r"carcinoma", r"\bhcc\b"
        ],
    ) else "NO"
)
m["Flag_HFD_signal"] = m["_combined"].apply(
    lambda x: "YES" if contains_any(
        x, [r"\bhfd\b", r"high[- ]fat", r"d12492", r"60%.*fat"]
    ) else "NO"
)
m["Flag_control_diet_signal"] = m["_combined"].apply(
    lambda x: "YES" if contains_any(
        x,
        [
            r"\bchow\b", r"\blfd\b", r"low[- ]fat", r"standard chow",
            r"\bscd\b", r"normal diet", r"control diet"
        ],
    ) else "NO"
)
m["Flag_special_diet_signal"] = m["_combined"].apply(
    lambda x: "YES" if contains_any(
        x,
        [
            r"cda[- ]?hfd", r"choline[- ]deficient", r"methionine.*choline",
            r"\bmcd\b", r"western diet", r"amln", r"\bgan\b",
            r"hfhc", r"ketogenic"
        ],
    ) else "NO"
)
m["Flag_intervention_signal"] = m["_combined"].apply(
    lambda x: "YES" if contains_any(
        x,
        [
            r"treated", r"treatment", r"drug", r"vehicle", r"compound",
            r"agonist", r"inhibitor", r"tirzepatide", r"exercise",
            r"reversal", r"switch"
        ],
    ) else "NO"
)
m["Flag_genetic_signal"] = m["_combined"].apply(
    lambda x: "YES" if contains_any(
        x,
        [
            r"knockout", r"\bko\b", r"transgenic", r"mutant",
            r"flox", r"cre\b", r"overexpress", r"deficient",
            r"null", r"heterozyg"
        ],
    ) else "NO"
)

def review_status(r):
    if r["Flag_mouse"] != "YES":
        return "REVIEW_EXCLUDE_ORGANISM"
    if r["Flag_bulk_RNAseq"] != "YES":
        return "REVIEW_EXCLUDE_NON_BULK_RNASEQ"
    if r["Raw_primary_status"] == "EXCLUDE_DATASET_FROM_RAW_PRIMARY":
        return "RAW_UNAVAILABLE_DATASET"
    if r["Raw_primary_status"] == "EXCLUDE_SAMPLE":
        return "RAW_UNAVAILABLE_SAMPLE"
    if r["Flag_liver_text"] != "YES":
        return "REVIEW_TISSUE"
    if r["Flag_non_whole_liver_signal"] == "YES":
        return "REVIEW_WHOLE_LIVER"
    if r["Flag_special_diet_signal"] == "YES":
        return "REVIEW_SPECIAL_DIET"
    if r["Flag_HFD_signal"] == "YES":
        return "HFD_CANDIDATE"
    if r["Flag_control_diet_signal"] == "YES":
        return "CONTROL_CANDIDATE"
    return "REVIEW_METADATA"

m["Sample_review_status_v1_1"] = m.apply(review_status, axis=1)
# legacy field replaced with corrected triage for convenience
m["Sample_review_status_v1"] = m["Sample_review_status_v1_1"]

for c in [
    "Manual_include_primary", "Manual_role", "Manual_exclusion_reason",
    "Sex_std", "Strain_std", "Age_std", "Diet_std", "Diet_fat_percent",
    "HFD_duration_weeks", "Genotype_std", "Treatment_std",
    "Cohort_id", "Contrast_id", "Reviewer_notes"
]:
    if c not in m.columns:
        m[c] = ""
    else:
        m[c] = ""

m["Acquisition_master_version"] = "v1.1"
m["Acquisition_build_time"] = NOW
m = m.drop(columns=["_combined"])

m.to_csv(MASTER, sep="\t", index=False)

run_master = runs.copy()
run_master["Acquisition_table_version"] = "v1.1"
run_master["Acquisition_build_time"] = NOW
run_master.to_csv(RUN_MASTER, sep="\t", index=False)

g = (
    m.groupby("GSE")
    .agg(
        N_GSM=("GSM", "nunique"),
        N_bulk_RNAseq=("Flag_bulk_RNAseq", lambda x: int((x == "YES").sum())),
        N_raw_mapped=("Raw_primary_status", lambda x: int((x == "RAW_MAPPED").sum())),
        N_HFD_candidate=("Sample_review_status_v1_1", lambda x: int((x == "HFD_CANDIDATE").sum())),
        N_control_candidate=("Sample_review_status_v1_1", lambda x: int((x == "CONTROL_CANDIDATE").sum())),
        N_review_other=(
            "Sample_review_status_v1_1",
            lambda x: int((~x.isin(["HFD_CANDIDATE", "CONTROL_CANDIDATE"])).sum())
        ),
        GSE_raw_harmonized_status=("GSE_raw_harmonized_status", uniq_join),
    )
    .reset_index()
)
g["Needs_manual_review"] = "YES"
g.to_csv(GSE_SUMMARY, sep="\t", index=False)

review = m[
    ~m["Sample_review_status_v1_1"].isin(["HFD_CANDIDATE", "CONTROL_CANDIDATE"])
].copy()
review.to_csv(REVIEW_QUEUE, sep="\t", index=False)

# ------------------------------------------------------------------
# Human-readable review pack
# ------------------------------------------------------------------
display_cols = [
    c for c in [
        "GSM", "Sample_title", "Source_name", "Tissue", "Diet",
        "Genotype", "Sex", "Age", "Strain", "Treatment",
        "Library_strategy", "Raw_availability", "N_SRR",
        "Sample_review_status_v1_1"
    ] if c in m.columns
]

with REVIEW_PACK.open("w") as f:
    for idx, (gse, sub) in enumerate(m.groupby("GSE", sort=True), start=1):
        f.write("#" * 100 + "\n")
        f.write(f"### {gse}   [{idx}/107]\n")
        f.write("#" * 100 + "\n")
        counts = sub["Sample_review_status_v1_1"].value_counts().to_dict()
        f.write(
            "Status counts: " +
            " | ".join(f"{k}={v}" for k, v in sorted(counts.items())) +
            "\n\n"
        )
        f.write(sub[display_cols].to_string(index=False))
        f.write("\n\n")

rows = []
for gse, sub in m.groupby("GSE", sort=True):
    c = sub["Sample_review_status_v1_1"].value_counts().to_dict()
    rows.append({
        "GSE": gse,
        "N_GSM": sub["GSM"].nunique(),
        "N_HFD_CANDIDATE": c.get("HFD_CANDIDATE", 0),
        "N_CONTROL_CANDIDATE": c.get("CONTROL_CANDIDATE", 0),
        "N_REVIEW_WHOLE_LIVER": c.get("REVIEW_WHOLE_LIVER", 0),
        "N_REVIEW_SPECIAL_DIET": c.get("REVIEW_SPECIAL_DIET", 0),
        "N_REVIEW_TISSUE": c.get("REVIEW_TISSUE", 0),
        "N_REVIEW_METADATA": c.get("REVIEW_METADATA", 0),
        "N_RAW_UNAVAILABLE_SAMPLE": c.get("RAW_UNAVAILABLE_SAMPLE", 0),
        "N_RAW_UNAVAILABLE_DATASET": c.get("RAW_UNAVAILABLE_DATASET", 0),
        "N_NON_BULK_RNASEQ": c.get("REVIEW_EXCLUDE_NON_BULK_RNASEQ", 0),
    })
pd.DataFrame(rows).to_csv(REVIEW_SUMMARY, sep="\t", index=False)

# ------------------------------------------------------------------
# QA
# ------------------------------------------------------------------
n_bulk = int((m["Flag_bulk_RNAseq"] == "YES").sum())
n_raw = int((m["Raw_primary_status"] == "RAW_MAPPED").sum())
n_lost = int((m["Raw_primary_status"] == "EXCLUDE_SAMPLE").sum())
n_dataset_unavailable = int(
    (m["Raw_primary_status"] == "EXCLUDE_DATASET_FROM_RAW_PRIMARY").sum()
)

assert len(m) == 3535
assert m["GSE"].nunique() == 107
assert n_bulk == 3368, f"Expected 3368 bulk RNA-seq GSM, got {n_bulk}"
assert n_lost == 1
assert n_dataset_unavailable == 12
assert n_raw + n_lost + n_dataset_unavailable == n_bulk

MAP_SUMMARY.write_text(
    "\n".join([
        "PRIMARY 107 SAMPLE/RUN MAPPING V1.1",
        "====================================",
        "Correction: include ssRNA-seq as bulk PolyA RNA-seq",
        f"Target GSE: {m['GSE'].nunique()}",
        f"Total GSM: {len(m)}",
        f"Bulk RNA-seq GSM: {n_bulk}",
        f"Mapped public-SRA bulk RNA-seq GSM: {n_raw}",
        f"Raw-lost GSM: {n_lost}",
        f"Raw-unavailable dataset GSM: {n_dataset_unavailable}",
        f"Run rows: {len(runs)}",
        f"Unique SRR: {runs['SRR'].nunique()}",
        f"ssRNA-seq GSM rescued: {len(need_rescue)}",
        f"ssRNA-seq rescue unresolved: {len(unresolved_df)}",
        "",
        f"Created: {NOW}",
    ]) + "\n"
)

CORRECTION.write_text(
    "\n".join([
        "DATA ACQUISITION PARSER CORRECTION V1.1",
        "=======================================",
        f"Date: {NOW}",
        "",
        "Issue:",
        'v1 classified bulk RNA-seq using exact Library_strategy == "RNA-Seq".',
        'GSE296500 uses Library_strategy == "ssRNA-seq" for 45 whole-liver PolyA RNA-seq samples.',
        "",
        "Correction:",
        '- Treat "RNA-Seq" and "ssRNA-seq" as eligible bulk RNA-seq strategies.',
        "- Rescue the 45 pre-existing SRX accessions to public SRA runs.",
        "- Rebuild acquisition/review outputs as v1.1.",
        "",
        "Scientific screening freeze affected: NO",
        "Primary GSE universe affected: NO (remains 107)",
        "Raw-phase mapping affected: YES",
        "Original v1 files modified: NO",
        "",
        "External verification:",
        "NCBI GEO GSE296500 describes the Series as PolyA RNA-seq and states that raw data are available in SRA.",
    ]) + "\n"
)

manifest = [
    "DISSERTATION-GRADE DATA ACQUISITION MANIFEST V1.1",
    "==================================================",
    f"Build time: {NOW}",
    "Frozen primary GSE universe: 107",
    f"Total GSM rows: {len(m)}",
    f"Bulk RNA-seq GSM rows: {n_bulk}",
    f"Bulk RNA-seq GSM with mapped public SRA runs: {n_raw}",
    f"Raw-lost individual GSM exceptions: {n_lost}",
    f"Raw-unavailable dataset GSM rows: {n_dataset_unavailable}",
    f"Unique SRR in run master: {runs['SRR'].nunique()}",
    "",
    "Correction from v1:",
    '- Added 45 GSE296500 samples labeled "ssRNA-seq" to the bulk RNA-seq universe.',
    "",
    "INPUT SHA256",
    f"{GEO}: {sha256(GEO)}",
    f"{RUNS_V1}: {sha256(RUNS_V1)}",
]
if EXC.exists():
    manifest.append(f"{EXC}: {sha256(EXC)}")
if DEC.exists():
    manifest.append(f"{DEC}: {sha256(DEC)}")

MANIFEST.write_text("\n".join(manifest) + "\n")
with MANIFEST.open("a") as f:
    f.write("\nOUTPUT SHA256\n")
    for p in [
        RUNS_V11, MASTER, RUN_MASTER, GSE_SUMMARY,
        REVIEW_QUEUE, REVIEW_PACK, REVIEW_SUMMARY,
        MAP_SUMMARY, CORRECTION
    ]:
        f.write(f"{p}: {sha256(p)}\n")

print()
print("================================================")
print("V1.1 CORRECTION COMPLETE")
print("================================================")
print("GSE:", m["GSE"].nunique())
print("Total GSM:", len(m))
print("Bulk RNA-seq GSM:", n_bulk)
print("RAW_MAPPED bulk RNA-seq GSM:", n_raw)
print("Raw exceptions:", n_lost + n_dataset_unavailable)
print("Unique SRR:", runs["SRR"].nunique())
print("ssRNA-seq rescued GSM:", len(need_rescue))
print("Unresolved:", len(unresolved_df))
print()
print("Review status counts:")
print(m["Sample_review_status_v1_1"].value_counts().to_string())
print()
print("Outputs:")
for p in [
    RUNS_V11, MASTER, GSE_SUMMARY, REVIEW_SUMMARY,
    MANIFEST, MAP_SUMMARY, CORRECTION
]:
    print(" ", p)
