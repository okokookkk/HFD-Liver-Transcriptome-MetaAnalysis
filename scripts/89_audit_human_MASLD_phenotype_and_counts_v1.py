#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np
import tarfile
import gzip
import re
import csv

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

AUDITDIR = ROOT / "00_metadata/human_validation/cohort_audit_v1"
SUPPDIR = ROOT / "00_metadata/human_validation/supplementary_v1"

SAMPLE_META = AUDITDIR / "human_MASLD_public_cohort_sample_metadata_v1.tsv"
SUPP_INDEX = AUDITDIR / "human_MASLD_public_cohort_supplementary_files_v1.tsv"

OUTDIR = ROOT / "00_metadata/human_validation/constructability_audit_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_PHENO = OUTDIR / "human_MASLD_candidate_phenotype_field_values_v1.tsv"
OUT_FILES = OUTDIR / "human_MASLD_count_file_structure_audit_v1.tsv"
OUT_HEADERS = OUTDIR / "human_MASLD_count_matrix_sample_headers_v1.tsv"
OUT_G135 = OUTDIR / "GSE135251_perGSM_file_audit_v1.tsv"
OUT_SUMMARY = OUTDIR / "human_MASLD_constructability_audit_summary_v1.txt"

for q in [SAMPLE_META, SUPP_INDEX]:
    if not q.exists():
        raise RuntimeError(f"Missing required input: {q}")

meta = pd.read_csv(SAMPLE_META, sep="\t", dtype=str).fillna("")
supp = pd.read_csv(SUPP_INDEX, sep="\t", dtype=str).fillna("")

EXPECTED = {
    "GSE126848": {"n": 57, "case": 31, "control": 26},
    "GSE130970": {"n": 78, "case": 72, "control": 6},
    "GSE135251": {"n": 216, "case": 206, "control": 10},
}

# ---------------------------------------------------------------------
# 1) Candidate phenotype fields: preserve GEO text and show distributions.
# ---------------------------------------------------------------------
keywords = [
    "disease", "diagnos", "condition", "group", "phenotype",
    "histolog", "nash", "nafl", "nafld", "masld", "control",
    "steatos", "fibrosis", "nas", "stage", "status", "type",
    "source_name", "characteristics", "title"
]
kw_re = re.compile("|".join(map(re.escape, keywords)), flags=re.I)

pheno_rows = []

for gse, g in meta.groupby("GSE", sort=False):
    for col in g.columns:
        if col in {"GSE", "GSM"}:
            continue
        vals = g[col].astype(str).str.strip()
        nonblank = vals[vals != ""]
        if nonblank.empty:
            continue

        uniq = sorted(nonblank.unique().tolist())
        # Candidate if column name is relevant OR values mention phenotype terms.
        values_blob = " | ".join(uniq[:100])
        candidate = bool(kw_re.search(col) or kw_re.search(values_blob))
        if not candidate:
            continue

        vc = nonblank.value_counts(dropna=False)
        for val, n in vc.items():
            pheno_rows.append({
                "GSE": gse,
                "Field": col,
                "Value": val,
                "N_samples": int(n),
                "N_unique_values_in_field": int(nonblank.nunique()),
                "Field_nonblank_samples": int(len(nonblank)),
            })

pheno = pd.DataFrame(pheno_rows)
if len(pheno):
    pheno = pheno.sort_values(
        ["GSE", "Field", "N_samples", "Value"],
        ascending=[True, True, False, True]
    )
pheno.to_csv(OUT_PHENO, sep="\t", index=False)

# ---------------------------------------------------------------------
# Helpers for count files.
# ---------------------------------------------------------------------
def find_download(gse, filename_substr):
    z = supp[(supp["GSE"] == gse) & supp["Filename"].str.contains(
        filename_substr, regex=False, na=False
    )]
    if len(z) != 1:
        raise RuntimeError(
            f"{gse}: expected one supplementary file containing "
            f"{filename_substr!r}, found {len(z)}"
        )
    return Path(z.iloc[0]["Downloaded_path"])

def open_text_auto(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")

def sniff_delim(path):
    with open_text_auto(path) as fh:
        lines = []
        for _ in range(5):
            line = fh.readline()
            if not line:
                break
            lines.append(line)
    sample = "".join(lines)
    if "\t" in sample:
        return "\t"
    if "," in sample:
        return ","
    try:
        return csv.Sniffer().sniff(sample).delimiter
    except Exception:
        return None

def inspect_matrix(path, gse, expected_n):
    delim = sniff_delim(path)
    if delim is None:
        raise RuntimeError(f"{gse}: could not determine delimiter for {path}")

    df = pd.read_csv(path, sep=delim, nrows=20, dtype=str)
    header = list(df.columns)
    if len(header) < 2:
        raise RuntimeError(f"{gse}: count matrix has <2 columns")

    feature_col = header[0]
    sample_cols = header[1:]

    # Read full first column only to get feature row count.
    full_feature = pd.read_csv(path, sep=delim, usecols=[0], dtype=str)
    n_features = len(full_feature)
    n_sample_cols = len(sample_cols)

    # Numeric/integer-likeness audit using first 20 rows.
    num = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    numeric_fraction = float(num.notna().mean().mean()) if num.size else np.nan
    arr = num.to_numpy(dtype=float)
    finite = np.isfinite(arr)
    if finite.any():
        integer_like_fraction = float(
            (np.abs(arr[finite] - np.round(arr[finite])) < 1e-8).mean()
        )
        nonnegative_fraction = float((arr[finite] >= 0).mean())
    else:
        integer_like_fraction = np.nan
        nonnegative_fraction = np.nan

    return {
        "GSE": gse,
        "File": str(path.relative_to(ROOT)),
        "Format": "single_matrix",
        "Delimiter": repr(delim),
        "Feature_column": feature_col,
        "N_features": n_features,
        "N_sample_columns": n_sample_cols,
        "Expected_samples": expected_n,
        "Sample_column_count_matches_expected": "YES" if n_sample_cols == expected_n else "NO",
        "First20_numeric_fraction": numeric_fraction,
        "First20_integer_like_fraction": integer_like_fraction,
        "First20_nonnegative_fraction": nonnegative_fraction,
    }, sample_cols

# ---------------------------------------------------------------------
# 2) GSE126848 and GSE130970 count matrix structures.
# ---------------------------------------------------------------------
file_rows = []
header_rows = []

g126 = find_download("GSE126848", "GSE126848_Gene_counts_raw.txt.gz")
r, hdr = inspect_matrix(g126, "GSE126848", EXPECTED["GSE126848"]["n"])
file_rows.append(r)
for i, s in enumerate(hdr, start=1):
    header_rows.append({
        "GSE": "GSE126848",
        "Column_index": i,
        "Count_matrix_sample_header": s,
        "Header_contains_GSM": re.search(r"GSM\d+", s).group(0)
            if re.search(r"GSM\d+", s) else "",
    })

g130 = find_download(
    "GSE130970",
    "GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz"
)
r, hdr = inspect_matrix(g130, "GSE130970", EXPECTED["GSE130970"]["n"])
file_rows.append(r)
for i, s in enumerate(hdr, start=1):
    header_rows.append({
        "GSE": "GSE130970",
        "Column_index": i,
        "Count_matrix_sample_header": s,
        "Header_contains_GSM": re.search(r"GSM\d+", s).group(0)
            if re.search(r"GSM\d+", s) else "",
    })

# ---------------------------------------------------------------------
# 3) GSE135251 RAW.tar contents and per-GSM count file audit.
# ---------------------------------------------------------------------
g135tar = find_download("GSE135251", "GSE135251_RAW.tar")

g135_rows = []
with tarfile.open(g135tar, "r") as tf:
    members = [m for m in tf.getmembers() if m.isfile()]

    for m in members:
        name = Path(m.name).name
        gsm_match = re.search(r"(GSM\d+)", name)
        gsm = gsm_match.group(1) if gsm_match else ""

        first_line = ""
        second_line = ""
        inspect_status = "NOT_TEXT_INSPECTED"

        try:
            fobj = tf.extractfile(m)
            raw = fobj.read(65536) if fobj is not None else b""

            if name.endswith(".gz"):
                try:
                    raw = gzip.decompress(raw)
                except Exception:
                    # gzip member may be larger than buffer; stream-decompress
                    fobj = tf.extractfile(m)
                    with gzip.GzipFile(fileobj=fobj, mode="rb") as gz:
                        raw = gz.read(65536)

            txt = raw.decode("utf-8", errors="replace")
            lines = txt.splitlines()
            if lines:
                first_line = lines[0][:500]
            if len(lines) > 1:
                second_line = lines[1][:500]
            inspect_status = "TEXT_PREVIEW_OK"
        except Exception as e:
            inspect_status = f"PREVIEW_FAIL:{type(e).__name__}"

        g135_rows.append({
            "Member_name": m.name,
            "Filename": name,
            "Size_bytes": int(m.size),
            "GSM_from_filename": gsm,
            "Preview_status": inspect_status,
            "First_line": first_line,
            "Second_line": second_line,
        })

g135 = pd.DataFrame(g135_rows)
g135.to_csv(OUT_G135, sep="\t", index=False)

n_g135_files = len(g135)
n_g135_gsm = g135.loc[g135["GSM_from_filename"] != "", "GSM_from_filename"].nunique()

file_rows.append({
    "GSE": "GSE135251",
    "File": str(g135tar.relative_to(ROOT)),
    "Format": "RAW_tar_per_sample_files",
    "Delimiter": "see_perGSM_audit",
    "Feature_column": "see_perGSM_audit",
    "N_features": np.nan,
    "N_sample_columns": n_g135_files,
    "Expected_samples": EXPECTED["GSE135251"]["n"],
    "Sample_column_count_matches_expected": "YES" if n_g135_gsm == EXPECTED["GSE135251"]["n"] else "NO",
    "First20_numeric_fraction": np.nan,
    "First20_integer_like_fraction": np.nan,
    "First20_nonnegative_fraction": np.nan,
})

files_audit = pd.DataFrame(file_rows)
files_audit.to_csv(OUT_FILES, sep="\t", index=False)

headers = pd.DataFrame(header_rows)
headers.to_csv(OUT_HEADERS, sep="\t", index=False)

# ---------------------------------------------------------------------
# 4) GSM/header mapping heuristics, WITHOUT assigning phenotype.
# ---------------------------------------------------------------------
summary_lines = [
    "===== HUMAN MASLD PHENOTYPE / COUNT CONSTRUCTABILITY AUDIT v1 =====",
    "",
    "STATUS=PASS_FOR_MANUAL_GROUP_ADJUDICATION",
    "HUMAN_DE_ANALYSIS_PERFORMED=NO",
    "CASE_CONTROL_ASSIGNMENT_PERFORMED=NO",
    "",
]

for gse in EXPECTED:
    gm = meta[meta["GSE"] == gse]
    summary_lines.append(
        f"{gse}_GEO_samples={gm['GSM'].nunique()};"
        f" expected={EXPECTED[gse]['n']}"
    )

summary_lines += [
    "",
    f"GSE126848_count_sample_columns={files_audit.loc[files_audit.GSE=='GSE126848','N_sample_columns'].iloc[0]}",
    f"GSE130970_count_sample_columns={files_audit.loc[files_audit.GSE=='GSE130970','N_sample_columns'].iloc[0]}",
    f"GSE135251_tar_files={n_g135_files}",
    f"GSE135251_unique_GSM_in_filenames={n_g135_gsm}",
    "",
    "EXPECTED_PRIMARY_CONTRASTS:",
    "GSE126848=31_NAFL_or_NASH_vs_26_histologically_normal",
    "GSE130970=72_NAFLD_spectrum_vs_6_histologically_normal",
    "GSE135251=206_NAFLD_vs_10_control",
    "",
    "MANUAL_GROUP_CONFIRMATION_REQUIRED=YES",
    "DO_NOT_INFER_GROUP_FROM_SAMPLE_TITLE_ALONE=YES",
    "",
    f"Candidate_phenotype_values={OUT_PHENO.relative_to(ROOT)}",
    f"Count_file_structure={OUT_FILES.relative_to(ROOT)}",
    f"Count_matrix_headers={OUT_HEADERS.relative_to(ROOT)}",
    f"GSE135251_perGSM_audit={OUT_G135.relative_to(ROOT)}",
]

OUT_SUMMARY.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
print("\n".join(summary_lines))

print("\n===== CANDIDATE PHENOTYPE FIELDS: COMPACT VIEW =====")
if pheno.empty:
    print("NONE")
else:
    # Print only fields with <=15 distinct values so output is readable.
    compact = pheno[pheno["N_unique_values_in_field"] <= 15].copy()
    for gse in EXPECTED:
        print(f"\n--- {gse} ---")
        z = compact[compact["GSE"] == gse]
        if z.empty:
            print("No compact candidate phenotype fields.")
            continue
        for field, q in z.groupby("Field", sort=False):
            vals = "; ".join(
                f"{r.Value} [{r.N_samples}]"
                for r in q.itertuples(index=False)
            )
            print(f"{field}: {vals}")

print("\n===== COUNT FILE STRUCTURE =====")
print(files_audit.to_string(index=False))

