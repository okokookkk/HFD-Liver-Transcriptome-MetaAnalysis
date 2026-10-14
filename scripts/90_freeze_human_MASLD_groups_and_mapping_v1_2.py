#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import re
import hashlib

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

META = ROOT / "00_metadata/human_validation/cohort_audit_v1/human_MASLD_public_cohort_sample_metadata_v1.tsv"
HEADERS = ROOT / "00_metadata/human_validation/constructability_audit_v1/human_MASLD_count_matrix_sample_headers_v1.tsv"

V1DIR = ROOT / "00_metadata/human_validation/frozen_groups_v1"
V1_GROUPS = V1DIR / "human_MASLD_primary_sample_groups_frozen_v1.tsv"
V1_G130_SENS = V1DIR / "GSE130970_strict4_control_sensitivity_groups_v1.tsv"
V1_MAP = V1DIR / "human_MASLD_count_sample_mapping_audit_v1.tsv"

OUTDIR = ROOT / "00_metadata/human_validation/frozen_groups_v1_2"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_GROUPS = OUTDIR / "human_MASLD_primary_sample_groups_frozen_v1_2.tsv"
OUT_G130_SENS = OUTDIR / "GSE130970_strict4_control_sensitivity_groups_v1_2.tsv"
OUT_MAP = OUTDIR / "human_MASLD_count_sample_mapping_audit_v1_2.tsv"
OUT_G126 = OUTDIR / "GSE126848_internal_id_crosswalk_v1_2.tsv"
OUT_SUM = OUTDIR / "human_MASLD_group_freeze_summary_v1_2.txt"

for q in [META, HEADERS, V1_GROUPS, V1_G130_SENS, V1_MAP]:
    if not q.exists():
        raise RuntimeError(f"Missing required input: {q}")

meta = pd.read_csv(META, sep="\t", dtype=str).fillna("")
headers = pd.read_csv(HEADERS, sep="\t", dtype=str).fillna("")
groups_v1 = pd.read_csv(V1_GROUPS, sep="\t", dtype=str).fillna("")
g130_sens_v1 = pd.read_csv(V1_G130_SENS, sep="\t", dtype=str).fillna("")
map_v1 = pd.read_csv(V1_MAP, sep="\t", dtype=str).fillna("")

def norm_text(x):
    return re.sub(r"[^A-Za-z0-9]+", "", str(x)).lower()

def numeric_key(x):
    s = str(x).strip()
    if re.fullmatch(r"\d+", s):
        return str(int(s))
    return None

def membership_sha(df):
    cols = ["GSE", "GSM", "Human_group"]
    z = df[cols].copy().sort_values(cols).astype(str)
    blob = "\n".join("\t".join(row) for row in z.itertuples(index=False, name=None))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()

groups_v1.to_csv(OUT_GROUPS, sep="\t", index=False)
g130_sens_v1.to_csv(OUT_G130_SENS, sep="\t", index=False)

sha_before = membership_sha(groups_v1)
sha_after = membership_sha(pd.read_csv(OUT_GROUPS, sep="\t", dtype=str).fillna(""))
if sha_before != sha_after:
    raise RuntimeError("Frozen biological membership changed")

g126m = meta[meta.GSE == "GSE126848"].copy()
g126h = headers[headers.GSE == "GSE126848"].copy()

desc_cols = [
    c for c in g126m.columns
    if c.lower() == "description" or c.lower().startswith("description.")
]
if not desc_cols:
    raise RuntimeError("No GSE126848 Description column found")

text_alias = {}
num_alias = {}
source_alias = {}

def add_text_alias(alias, gsm, source):
    k = norm_text(alias)
    if not k:
        return
    text_alias.setdefault(k, set()).add(gsm)
    source_alias.setdefault(("TEXT", k, gsm), set()).add(source)

def add_num_alias(alias, gsm, source):
    k = numeric_key(alias)
    if k is None:
        return
    num_alias.setdefault(k, set()).add(gsm)
    source_alias.setdefault(("NUM", k, gsm), set()).add(source)

for _, r in g126m.iterrows():
    gsm = r["GSM"]
    add_text_alias(gsm, gsm, "GSM")
    if "title" in g126m.columns:
        add_text_alias(r["title"], gsm, "GEO_title")
    for c in desc_cols:
        add_text_alias(r[c], gsm, f"GEO_{c}_exact")
        add_num_alias(r[c], gsm, f"GEO_{c}_numeric_equivalence")

rows = []
for _, r in g126h.iterrows():
    header = str(r["Count_matrix_sample_header"]).strip()
    exact_candidates = sorted(text_alias.get(norm_text(header), set()))
    numk = numeric_key(header)
    numeric_candidates = sorted(num_alias.get(numk, set())) if numk is not None else []

    if len(exact_candidates) == 1:
        gsm = exact_candidates[0]
        status = "PASS_EXACT_METADATA_ALIAS"
        method = ";".join(sorted(source_alias.get(("TEXT", norm_text(header), gsm), set())))
    elif len(exact_candidates) > 1:
        gsm = ""
        status = "AMBIGUOUS_EXACT"
        method = ""
    elif len(numeric_candidates) == 1:
        gsm = numeric_candidates[0]
        status = "PASS_NUMERIC_INTERNAL_ID_EQUIVALENCE"
        method = ";".join(sorted(source_alias.get(("NUM", numk, gsm), set())))
    elif len(numeric_candidates) > 1:
        gsm = ""
        status = "AMBIGUOUS_NUMERIC"
        method = ""
    else:
        m = re.search(r"(GSM\d+)", header)
        if m and m.group(1) in set(g126m.GSM):
            gsm = m.group(1)
            status = "PASS_EMBEDDED_GSM"
            method = "embedded_GSM"
        else:
            gsm = ""
            status = "UNMAPPED"
            method = ""

    if gsm:
        rr = g126m[g126m.GSM == gsm].iloc[0]
        desc = rr[desc_cols[0]]
        title = rr["title"] if "title" in rr.index else ""
        disease = rr["disease:ch1"] if "disease:ch1" in rr.index else (
            rr["characteristics_ch1.2"] if "characteristics_ch1.2" in rr.index else ""
        )
    else:
        desc = title = disease = ""

    rows.append({
        "GSE": "GSE126848",
        "Count_sample_identifier": header,
        "Mapped_GSM": gsm,
        "Mapping_status": status,
        "Mapping_method": method,
        "GEO_Description": desc,
        "GEO_title": title,
        "GEO_disease": disease,
        "Header_numeric_key": numeric_key(header) or "",
        "Description_numeric_key": numeric_key(desc) or "",
    })

g126_map = pd.DataFrame(rows)
g126_map.to_csv(OUT_G126, sep="\t", index=False)

other = map_v1[map_v1.GSE.isin(["GSE130970", "GSE135251"])].copy()
if "Mapping_method" not in other.columns:
    other["Mapping_method"] = other["Mapping_status"].map({
        "PASS_EXACT_ALIAS": "v1_exact_alias",
        "PASS_FILENAME_GSM": "filename_GSM",
    }).fillna("v1_preserved")

all_cols = list(dict.fromkeys(list(g126_map.columns) + list(other.columns)))
for c in all_cols:
    if c not in g126_map.columns:
        g126_map[c] = ""
    if c not in other.columns:
        other[c] = ""

mapping = pd.concat([g126_map[all_cols], other[all_cols]], ignore_index=True)
mapping.to_csv(OUT_MAP, sep="\t", index=False)

expected = {"GSE126848":57, "GSE130970":78, "GSE135251":216}
audit_rows = []
full_mapping = True

for gse, n_expected in expected.items():
    z = mapping[mapping.GSE == gse].copy()
    mapped = z.Mapped_GSM.astype(str) != ""
    n_count = len(z)
    n_mapped = int(mapped.sum())
    n_unique = int(z.loc[mapped, "Mapped_GSM"].nunique())
    n_unmapped = int((z.Mapping_status == "UNMAPPED").sum())
    n_ambiguous = int(z.Mapping_status.astype(str).str.startswith("AMBIGUOUS").sum())

    frozen = set(groups_v1.loc[groups_v1.GSE == gse, "GSM"])
    observed = set(z.loc[mapped, "Mapped_GSM"])
    exact_membership = frozen == observed

    ok = (
        n_count == n_expected and
        n_mapped == n_expected and
        n_unique == n_expected and
        n_unmapped == 0 and
        n_ambiguous == 0 and
        exact_membership
    )
    full_mapping = full_mapping and ok

    audit_rows.append({
        "GSE": gse,
        "N_count": n_count,
        "N_mapped": n_mapped,
        "N_unique": n_unique,
        "N_unmapped": n_unmapped,
        "N_ambiguous": n_ambiguous,
        "Frozen_membership_exact_match": exact_membership,
        "PASS": ok,
    })

audit = pd.DataFrame(audit_rows)

desc_numeric = []
for c in desc_cols:
    desc_numeric += [numeric_key(v) for v in g126m[c]]
desc_numeric = {v for v in desc_numeric if v is not None}
header_numeric = {numeric_key(v) for v in g126h.Count_matrix_sample_header}
header_numeric.discard(None)
numeric_set_match = desc_numeric == header_numeric

n_exact = int((g126_map.Mapping_status == "PASS_EXACT_METADATA_ALIAS").sum())
n_numeric = int((g126_map.Mapping_status == "PASS_NUMERIC_INTERNAL_ID_EQUIVALENCE").sum())
n_embedded = int((g126_map.Mapping_status == "PASS_EMBEDDED_GSM").sum())

group_counts = groups_v1.groupby(["GSE","Human_group"]).size().unstack(fill_value=0)

status = "PASS" if full_mapping and numeric_set_match else "FAIL"

summary = [
    "===== HUMAN MASLD GROUP FREEZE / COUNT MAPPING SUMMARY v1.2 =====",
    "",
    f"STATUS={status}",
    "SUPERSEDES=v1_and_v1_1_MAPPING_ONLY",
    "BIOLOGICAL_GROUP_MEMBERSHIP_CHANGED=NO",
    f"PRIMARY_GROUP_MEMBERSHIP_SHA256={sha_after}",
    "",
    "PRIMARY_HUMAN_EXTERNAL_VALIDATION_GROUPS:",
    f"GSE126848={int(group_counts.loc['GSE126848','CASE'])}_case_vs_{int(group_counts.loc['GSE126848','CONTROL'])}_control",
    f"GSE130970={int(group_counts.loc['GSE130970','CASE'])}_case_vs_{int(group_counts.loc['GSE130970','CONTROL'])}_control",
    f"GSE135251={int(group_counts.loc['GSE135251','CASE'])}_case_vs_{int(group_counts.loc['GSE135251','CONTROL'])}_control",
    "",
    "GSE130970_PRIMARY_CONTROL_N=6",
    "GSE130970_STRICT_SENSITIVITY_CONTROL_N=4",
    "",
    "GSE126848_MAPPING_REPAIR=GEO_Description_internal_sample_ID_with_numeric_leading_zero_equivalence",
    "GSE126848_MAPPING_BY_COLUMN_ORDER=NO",
    f"GSE126848_description_columns={';'.join(desc_cols)}",
    f"GSE126848_numeric_internal_ID_set_exact_match={'YES' if numeric_set_match else 'NO'}",
    f"GSE126848_exact_alias_mappings={n_exact}",
    f"GSE126848_numeric_equivalence_mappings={n_numeric}",
    f"GSE126848_embedded_GSM_mappings={n_embedded}",
    "",
    "COUNT_SAMPLE_MAPPING:",
]
for _, r in audit.iterrows():
    summary.append(
        f"{r.GSE}=count_samples:{int(r.N_count)};"
        f"mapped:{int(r.N_mapped)};"
        f"unique_GSM:{int(r.N_unique)};"
        f"unmapped:{int(r.N_unmapped)};"
        f"ambiguous:{int(r.N_ambiguous)};"
        f"frozen_membership_match:{'YES' if r.Frozen_membership_exact_match else 'NO'}"
    )

summary += [
    "",
    "HUMAN_DE_ANALYSIS_PERFORMED=NO",
    "GROUP_DEFINITIONS_FROZEN_BEFORE_HUMAN_DE=YES",
    f"COUNT_MEMBERSHIP_FROZEN={'YES' if status == 'PASS' else 'NO'}",
    f"Primary_groups={OUT_GROUPS.relative_to(ROOT)}",
    f"GSE130970_strict4_sensitivity={OUT_G130_SENS.relative_to(ROOT)}",
    f"Count_mapping_audit={OUT_MAP.relative_to(ROOT)}",
    f"GSE126848_crosswalk={OUT_G126.relative_to(ROOT)}",
]

OUT_SUM.write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))

if status != "PASS":
    raise SystemExit(2)

