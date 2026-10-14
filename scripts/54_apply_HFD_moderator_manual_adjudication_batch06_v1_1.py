#!/usr/bin/env python3
from pathlib import Path
import argparse, re
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=__import__("os").environ.get("HFD_PROJECT_ROOT", ".")
)
args = ap.parse_args()

ROOT = Path(args.root)
OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"

MASTER = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch05_v1.tsv"
GSM_AUTO = (
    ROOT /
    "00_metadata/contrasts/moderators/"
    "HFD_primary837_GSM_moderator_extraction_GEOauto_v1_1_1.tsv"
)
DURATION_POLICY = OUTDIR / "HFD_duration_unit_standardization_policy_v1.txt"

OUT = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch06_v1_1.tsv"
ATTEMPT_OUT = OUTDIR / "HFD_batch06_source_insufficient_review_log_v1_1.tsv"

for p in [MASTER, GSM_AUTO, DURATION_POLICY]:
    if not p.exists():
        raise FileNotFoundError(p)

x = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
g = pd.read_csv(GSM_AUTO, sep="\t", dtype=str).fillna("")

assert x["GSE"].nunique() == 64
assert x["Contrast_id"].nunique() == 105
assert g["GSM"].nunique() == 837

def fill_var(
    idx, var, manual_value, source_type, citation, evidence,
    confidence="HIGH", status="MANUAL_CONFIRMED", final_note=""
):
    x.loc[idx, f"{var}_manual_value"] = str(manual_value)
    x.loc[idx, f"{var}_manual_source_type"] = source_type
    x.loc[idx, f"{var}_manual_source_citation"] = citation
    x.loc[idx, f"{var}_manual_evidence"] = evidence
    x.loc[idx, f"{var}_manual_confidence"] = confidence
    x.loc[idx, f"{var}_adjudication_status"] = status
    x.loc[idx, f"{var}_final"] = str(manual_value)
    x.loc[idx, f"{var}_final_confidence"] = confidence
    x.loc[idx, f"{var}_final_source"] = source_type
    x.loc[idx, f"{var}_final_note"] = final_note

def mark_not_reported(idx, var, source_type, citation, evidence):
    x.loc[idx, f"{var}_manual_value"] = ""
    x.loc[idx, f"{var}_manual_source_type"] = source_type
    x.loc[idx, f"{var}_manual_source_citation"] = citation
    x.loc[idx, f"{var}_manual_evidence"] = evidence
    x.loc[idx, f"{var}_manual_confidence"] = "UNRESOLVED"
    x.loc[idx, f"{var}_adjudication_status"] = "NOT_REPORTED"
    x.loc[idx, f"{var}_final"] = ""
    x.loc[idx, f"{var}_final_confidence"] = ""
    x.loc[idx, f"{var}_final_source"] = ""
    x.loc[idx, f"{var}_final_note"] = evidence

def set_pub(gse, pmid, doi, title):
    idx = x["GSE"].eq(gse)
    x.loc[idx, "publication_pmid"] = pmid
    x.loc[idx, "publication_doi"] = doi
    x.loc[idx, "publication_title"] = title
    return idx

def selected_text(gse):
    z = g[g["GSE"].eq(gse)].copy()
    vals = []
    for c in ["sample_title_raw", "sample_characteristics_raw"]:
        if c in z.columns:
            vals.extend(v for v in z[c].tolist() if v)
    return " || ".join(vals)

attempts = []

# ============================================================
# GSE293219 — COMPLETE, strict duration upgrade
# Tokgozoglu et al. Cell Death Dis. 2025
# PMID 41053093; DOI 10.1038/s41419-025-08004-1
#
# - male C57BL/6J
# - 16 weeks
# - HFD D12492, 60% energy from fat
# - standard chow NCD
# Frozen primary must be WT liver NCD vs WT liver HFD, not MlklAdi-KO.
# ============================================================
gse = "GSE293219"
idx = set_pub(
    gse,
    "41053093",
    "10.1038/s41419-025-08004-1",
    "Adipocyte-specific Mlkl knockout mitigates obesity-induced metabolic dysfunction by enhancing mitochondrial functions",
)
citation = (
    "Tokgozoglu et al. Cell Death Dis. 2025;16:683. "
    "PMID:41053093; DOI:10.1038/s41419-025-08004-1"
)
t = selected_text(gse)

if not re.search(r"\bLV\b|liver", t, re.I):
    raise RuntimeError("GSE293219: liver samples not detected.")
if not re.search(r"\bWT\b", t, re.I):
    raise RuntimeError("GSE293219: WT samples not detected.")
if not re.search(r"\bhfd\b", t, re.I):
    raise RuntimeError("GSE293219: HFD samples not detected.")
if not re.search(r"\bncd\b|chow", t, re.I):
    raise RuntimeError("GSE293219: NCD/chow samples not detected.")
if re.search(r"MlklAdi[-_ ]?KO|Adi[-_ ]?KO", t, re.I):
    raise RuntimeError(
        "GSE293219: Mlkl adipocyte-KO samples detected in frozen primary contrast."
    )

fill_var(
    idx, "sex", "MALE", "PUBLICATION_METHODS", citation,
    "Sixteen-week-old WT and MlklAdi-KO male mice were used for diet studies."
)
fill_var(
    idx, "strain", "C57BL/6J", "PUBLICATION_METHODS", citation,
    "Mlkl conditional mice and WT littermates were bred on a C57BL/6J genetic background."
)
fill_var(
    idx, "hfd_duration_weeks", "16", "PUBLICATION_METHODS", citation,
    "WT and MlklAdi-KO mice were fed NCD or HFD for 16 weeks."
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60", "PUBLICATION_METHODS", citation,
    "HFD was D12492 and provided 60% of total energy from lipid."
)
fill_var(
    idx, "control_diet_type", "CHOW", "PUBLICATION_METHODS", citation,
    "Control diet was standard chow NCD (#Rod 16-R, Genobios)."
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected samples verified as WT liver NCD vs WT liver HFD; "
    "MlklAdi-KO samples are rejected by QA."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH06_ADJUDICATED"

# ============================================================
# GSE268200 — COMPLETE, strict duration upgrade + duration correction
# Lopez-Pascual et al. Clin Sci (Lond). 2024
# PMID 39301694; DOI 10.1042/CS20241137
#
# - male C57BL/6NTac
# - DIO-MASH/AMLN: 40 kcal% fat, 20 kcal% fructose, 2% cholesterol
# - paired-age normal chow control
#
# Critical duration definition:
# publication reports 21 weeks HFD BEFORE the 6-week mRNA-LNP phase,
# with HFD continued during treatment; Methods explicitly identify
# terminal collection as 28 weeks on HFD.
#
# Frozen primary is the control-mRNA/vehicle HFD-vs-chow baseline,
# not NM-FGF19 intervention.
# ============================================================
gse = "GSE268200"
idx = set_pub(
    gse,
    "39301694",
    "10.1042/CS20241137",
    "Non-mitogenic FGF19 mRNA-based therapy for the treatment of experimental metabolic dysfunction-associated steatotic liver disease (MASLD)",
)
citation = (
    "Lopez-Pascual et al. Clin Sci (Lond). 2024;138:1265-1284. "
    "PMID:39301694; DOI:10.1042/CS20241137"
)
t = selected_text(gse)
z268 = g[g["GSE"].eq(gse)].copy()

# v1.1 QA correction:
# The frozen GSE268200 sample-specific title/characteristics do not reliably
# contain literal HFD/DIO/MASH or chow/lean tokens. Do not re-derive the
# already-frozen biological group assignment from those free-text strings.
#
# Instead verify the frozen sample universe itself (8 GSM, one contrast)
# and guard only against explicit NM-FGF19 intervention contamination.
if z268["GSM"].nunique() != 8:
    raise RuntimeError(
        f"GSE268200: expected 8 frozen GSM, observed {z268['GSM'].nunique()}."
    )
if z268["Contrast_id"].nunique() != 1:
    raise RuntimeError(
        f"GSE268200: expected one frozen contrast, observed "
        f"{z268['Contrast_id'].nunique()}."
    )
if not t.strip():
    raise RuntimeError("GSE268200: selected sample-specific metadata text is empty.")
if re.search(r"NM[-_ ]?FGF19|FGF19[-_ ]?mRNA", t, re.I):
    raise RuntimeError(
        "GSE268200: NM-FGF19 intervention samples detected in frozen primary contrast."
    )

print("\n===== GSE268200 FROZEN SAMPLE QA =====")
print("selected_GSM:", z268["GSM"].nunique())
print("selected_contrasts:", z268["Contrast_id"].nunique())
print("sample_specific_text_preview:", t[:1200])

fill_var(
    idx, "sex", "MALE", "PUBLICATION_METHODS", citation,
    "The DIO-MASH model used C57BL/6NTac male mice."
)
fill_var(
    idx, "strain", "C57BL/6NTac", "PUBLICATION_METHODS", citation,
    "The DIO-MASH model is explicitly reported as C57BL/6NTac."
)
fill_var(
    idx, "hfd_duration_weeks", "27", "PUBLICATION_METHODS", citation,
    "The publication states that DIO-MASH mice were fed the diet for 21 weeks "
    "before control/NM-FGF19 mRNA-LNP treatment, followed by another six weeks "
    "during which diet feeding continued; total HFD exposure to collection is "
    "therefore 27 weeks.",
    status="PAPER_GEO_CONFLICT_RESOLVED",
    final_note=(
        "GEO-auto 21 weeks captured the pre-treatment HFD phase, not total "
        "HFD exposure to tissue collection. Under the prespecified definition "
        "(HFD initiation to tissue collection), final duration is 21 + 6 = 27 weeks. "
        "This supersedes the erroneous 28-week value in Batch06 v1."
    ),
)
fill_var(
    idx, "hfd_fat_kcal_percent", "40", "PUBLICATION_METHODS", citation,
    "DIO-MASH/AMLN diet contained 40 kcal% fat, 20 kcal% fructose, and 2% cholesterol."
)
fill_var(
    idx, "control_diet_type", "CHOW", "PUBLICATION_METHODS", citation,
    "Paired-age lean control mice were maintained on normal chow diet."
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen contrast is the non-FGF19 baseline HFD-vs-chow comparison; "
    "NM-FGF19 intervention samples are rejected by QA."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH06_ADJUDICATED"
x.loc[idx, "general_note"] = (
    "Special AMLN/DIO-MASH HFD (fat+fructose+cholesterol). Retain diet "
    "composition as biological context even though it remains eligible."
)

# ============================================================
# GSE191004 — PARTIAL, strict duration upgrade
# Anagnostopoulos et al. Cell Death Dis. 2022
# PMID 35436993; DOI 10.1038/s41419-022-04834-5
#
# - adult male
# - 4 weeks
# - SAFE 260 HF = 60% energy from fat
# - regular chow SAFE A04
# - Gabrg2tm1Wul/J +/+ vs F77I/F77I context
#
# Exact background/substrain of the selected RNA-seq +/+ cohort is not
# explicitly stated in the study Methods; do not infer it.
# Frozen primary must be +/+ chow vs +/+ HFD, not F77I/F77I.
# ============================================================
gse = "GSE191004"
idx = set_pub(
    gse,
    "35436993",
    "10.1038/s41419-022-04834-5",
    "An obesogenic feedforward loop involving PPARγ, acyl-CoA binding protein and GABAA receptor",
)
citation = (
    "Anagnostopoulos et al. Cell Death Dis. 2022;13:356. "
    "PMID:35436993; DOI:10.1038/s41419-022-04834-5"
)
t = selected_text(gse)

if not re.search(r"High[-_ ]?fat|HFD", t, re.I):
    raise RuntimeError("GSE191004: HFD samples not detected.")
if not re.search(r"Regular[-_ ]?chow|RCD|chow", t, re.I):
    raise RuntimeError("GSE191004: regular-chow control not detected.")
if re.search(r"F77I/F77I|F77I_F77I|homozygous.*F77I", t, re.I):
    raise RuntimeError(
        "GSE191004: Gabrg2 F77I/F77I mutant samples detected in frozen primary contrast."
    )

fill_var(
    idx, "sex", "MALE", "GEO_SERIES", citation + "; GSE191004",
    "GSE191004 describes adult male mice; study Methods also use 8-12-week-old male Gabrg2tm1Wul/J mice."
)
fill_var(
    idx, "hfd_duration_weeks", "4", "GEO_SERIES", citation + "; GSE191004",
    "GSE191004 explicitly states a 4-week regular-chow/HFD RNA-seq regime."
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60", "DIET_MANUFACTURER_CODE",
    citation + "; SAFE 260 HF U8978 Version 19",
    "Study Methods specify SAFE #260 HF; SAFE manufacturer documentation "
    "identifies 260 HF as 60% of energy from fat."
)
fill_var(
    idx, "control_diet_type", "CHOW", "PUBLICATION_METHODS", citation,
    "Study Methods specify regular SAFE #A04 chow as control."
)
mark_not_reported(
    idx, "strain", "PUBLICATION_METHODS", citation,
    "The study identifies the RNA-seq mice as Gabrg2tm1Wul/J +/+ or F77I/F77I "
    "but does not explicitly state the selected +/+ cohort's C57BL/6 J/N "
    "background/substrain. The allele stock name is retained as genotype "
    "context and is not converted into a strain label."
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected contrast is regular-chow vs HFD within the non-F77I "
    "control genotype; mutant samples are rejected by QA."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH06_PARTIAL_REVIEWED"

attempts.append({
    "GSE": "GSE191004",
    "unresolved_fields": "strain",
    "reason": (
        "Exact selected RNA-seq control-cohort background/substrain is not "
        "explicitly reported; no J/N inference from the Gabrg2 allele name."
    ),
})

# ============================================================
# GSE205846 — COMPLETE, strict duration upgrade
# Xia et al. Nat Commun. 2023
# PMID 37914694; DOI 10.1038/s41467-023-42785-3
#
# Selected frozen primary must be WT chow vs WT HFD, not ERRa KO.
# - male
# - C57BL/6N
# - 15 weeks
# - HFD TD.06414 = 60.3 kcal% fat
# - control TD.08806 = 10.5 kcal% fat
#
# Although publication calls TD.08806 "control chow", it is a defined
# custom low-fat control; final analytical category = LFD.
# ============================================================
gse = "GSE205846"
idx = set_pub(
    gse,
    "37914694",
    "10.1038/s41467-023-42785-3",
    "Hepatocyte FBXW7-dependent activity of nutrient-sensing nuclear receptors controls systemic energy homeostasis and NASH progression in male mice",
)
citation = (
    "Xia et al. Nat Commun. 2023;14:6982. "
    "PMID:37914694; DOI:10.1038/s41467-023-42785-3"
)
t = selected_text(gse)

if not re.search(r"\bWT\b", t, re.I):
    raise RuntimeError("GSE205846: WT samples not detected.")
if not re.search(r"high[- ]fat|HFD", t, re.I):
    raise RuntimeError("GSE205846: HFD samples not detected.")
if not re.search(r"chow|control", t, re.I):
    raise RuntimeError("GSE205846: control samples not detected.")
if re.search(r"ERR.?a.?KO|ERRα.?KO|Esrra.?KO", t, re.I):
    raise RuntimeError(
        "GSE205846: ERRα-KO samples detected in frozen primary contrast."
    )

fill_var(
    idx, "sex", "MALE", "PUBLICATION_METHODS", citation,
    "Mouse experiments and publication title identify the HFD RNA-seq cohort as male."
)
fill_var(
    idx, "strain", "C57BL/6N", "GEO_SERIES", citation + "; GSE205846",
    "GSE205846 explicitly states that WT and ERRα-KO HFD RNA-seq mice were on a C57BL/6N background."
)
fill_var(
    idx, "hfd_duration_weeks", "15", "PUBLICATION_METHODS", citation,
    "WT and ERRα-KO mice were fed HFD or control diet for 15 weeks starting at 6 weeks of age."
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60.3", "PUBLICATION_METHODS", citation,
    "HFD was Envigo/Teklad TD.06414, 60.3 kcal% fat."
)
fill_var(
    idx, "control_diet_type", "LFD", "PUBLICATION_METHODS", citation,
    "Control was Envigo/Teklad TD.08806, 10.5 kcal% fat.",
    status="MANUAL_CONFIRMED",
    final_note=(
        "The paper calls TD.08806 a control chow diet, but it is a defined "
        "custom 10.5 kcal% fat control. For cross-study analytical "
        "classification it is treated as LFD/purified control, not standard chow."
    ),
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected samples verified as WT control-diet vs WT HFD; ERRα-KO samples are rejected by QA."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH06_ADJUDICATED"

# ============================================================
# GSE285614 — COMPLETE metadata, but duration remains MEDIUM
# Guan et al. J Clin Invest. 2025
# PMID 41178716; DOI 10.1172/JCI190635
#
# Critical parser correction:
# auto duration=8 weeks captured AGE AT HFD START, not exposure duration.
# Bulk liver RNA-seq WT chow/HFD cohort was fed diets for 3 months.
#
# Under existing duration-unit policy:
# 3 calendar months -> 13.0446 weeks, MEDIUM confidence.
# Therefore this correction improves validity but does NOT make duration strict-ready.
# ============================================================
gse = "GSE285614"
idx = set_pub(
    gse,
    "41178716",
    "10.1172/JCI190635",
    "Adipocyte death promotes hepatic infiltration of S100A8+ macrophages and steatotic liver disease progression in mice",
)
citation = (
    "Guan et al. J Clin Invest. 2025;135:e190635. "
    "PMID:41178716; DOI:10.1172/JCI190635"
)
t = selected_text(gse)

if not re.search(r"Chow_WT|chow.*WT", t, re.I):
    raise RuntimeError("GSE285614: WT chow samples not detected.")
if not re.search(r"HFD_WT|HFD.*WT", t, re.I):
    raise RuntimeError("GSE285614: WT HFD samples not detected.")
if re.search(r"Bcl2TG|S100a8|Cx3cr1Cre", t, re.I):
    raise RuntimeError(
        "GSE285614: transgenic/conditional-KO samples detected in frozen primary contrast."
    )

MONTH_TO_WEEK = 365.25 / 12.0 / 7.0
duration_3mo = f"{3.0 * MONTH_TO_WEEK:.4f}"

fill_var(
    idx, "sex", "MALE", "PUBLICATION_RESULTS", citation,
    "Bulk liver RNA-seq chow/HFD WT and Bcl2AdTG cohorts were male."
)
fill_var(
    idx, "strain", "C57BL/6J", "PUBLICATION_METHODS", citation,
    "C57BL/6J mice were used; selected frozen primary samples are WT."
)
fill_var(
    idx, "hfd_duration_weeks", duration_3mo, "PUBLICATION_RESULTS", citation,
    "Bulk liver RNA-seq WT/Bcl2AdTG chow-vs-HFD cohort was fed diets for 3 months.",
    confidence="MEDIUM",
    status="AUTO_EXTRACTION_ERROR_CORRECTED",
    final_note=(
        f"GEO-auto 8 weeks captured age at HFD initiation, not exposure duration. "
        f"Source duration='3 months'; standardized={duration_3mo} weeks using "
        "HFD_duration_unit_standardization_policy_v1.txt. "
        "MEDIUM confidence; excluded from strict-duration sensitivity."
    ),
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60", "PUBLICATION_METHODS", citation,
    "HFD was Research Diets D12492, 60 kcal% fat."
)
fill_var(
    idx, "control_diet_type", "CHOW", "PUBLICATION_METHODS", citation,
    "Control was chow diet, reported as 10% kcal fat."
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected samples verified as Chow_WT vs HFD_WT; Bcl2TG and "
    "S100a8/Cx3cr1Cre genetic cohorts are rejected by QA."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH06_ADJUDICATED"

# ------------------------------------------------------------
# QA
# ------------------------------------------------------------
batch_gses = [
    "GSE293219", "GSE268200", "GSE191004",
    "GSE205846", "GSE285614",
]
b = x[x["GSE"].isin(batch_gses)].copy()

assert b["GSE"].nunique() == 5
assert b["Contrast_id"].nunique() == 5

# Complete all-moderator GSE
for gse in ["GSE293219", "GSE268200", "GSE205846", "GSE285614"]:
    z = b[b["GSE"].eq(gse)]
    for var in [
        "sex", "strain", "hfd_duration_weeks",
        "hfd_fat_kcal_percent", "control_diet_type"
    ]:
        assert (z[f"{var}_final"] != "").all(), (gse, var)
    assert set(z["overall_adjudication_status"]) == {"BATCH06_ADJUDICATED"}

# Partial strain-only unresolved
z = b[b["GSE"].eq("GSE191004")]
assert set(z["strain_final"]) == {""}
for var in ["sex", "hfd_duration_weeks", "hfd_fat_kcal_percent", "control_diet_type"]:
    assert (z[f"{var}_final"] != "").all(), var
assert set(z["overall_adjudication_status"]) == {"BATCH06_PARTIAL_REVIEWED"}

# Exact values
assert set(b[b["GSE"].eq("GSE293219")]["hfd_duration_weeks_final"]) == {"16"}
assert set(b[b["GSE"].eq("GSE293219")]["hfd_fat_kcal_percent_final"]) == {"60"}

assert set(b[b["GSE"].eq("GSE268200")]["hfd_duration_weeks_final"]) == {"27"}
assert set(b[b["GSE"].eq("GSE268200")]["hfd_fat_kcal_percent_final"]) == {"40"}
assert set(b[b["GSE"].eq("GSE268200")]["strain_final"]) == {"C57BL/6NTac"}

assert set(b[b["GSE"].eq("GSE191004")]["hfd_duration_weeks_final"]) == {"4"}
assert set(b[b["GSE"].eq("GSE191004")]["hfd_fat_kcal_percent_final"]) == {"60"}

assert set(b[b["GSE"].eq("GSE205846")]["hfd_duration_weeks_final"]) == {"15"}
assert set(b[b["GSE"].eq("GSE205846")]["hfd_fat_kcal_percent_final"]) == {"60.3"}
assert set(b[b["GSE"].eq("GSE205846")]["control_diet_type_final"]) == {"LFD"}

assert set(b[b["GSE"].eq("GSE285614")]["hfd_duration_weeks_final"]) == {duration_3mo}
assert set(b[b["GSE"].eq("GSE285614")]["hfd_duration_weeks_final_confidence"]) == {"MEDIUM"}

# Preserve previously reviewed statuses
prev = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
old_reviewed = prev[
    prev["overall_adjudication_status"].str.contains(
        "ADJUDICATED|PARTIAL_REVIEWED", regex=True, na=False
    )
]["GSE"].unique()

for gse in old_reviewed:
    a = prev[prev["GSE"].eq(gse)].sort_values("Contrast_id")
    bb = x[x["GSE"].eq(gse)].sort_values("Contrast_id")
    assert (
        a["overall_adjudication_status"].tolist()
        == bb["overall_adjudication_status"].tolist()
    ), gse

x.to_csv(OUT, sep="\t", index=False)
pd.DataFrame(attempts).to_csv(ATTEMPT_OUT, sep="\t", index=False)

CORRECTION_NOTE = OUTDIR / "HFD_batch06_v1_to_v1_1_correction_note.txt"
CORRECTION_NOTE.write_text(
    "Batch06 v1 -> v1.1 correction note\n"
    "\n"
    "GSE268200 corrections:\n"
    "1. QA guard correction: sample-specific GEO title/characteristics do not "
    "reliably contain literal HFD/DIO/MASH or chow/lean tokens. v1.1 therefore "
    "verifies the frozen sample universe (8 GSM, one contrast) and rejects only "
    "explicit NM-FGF19 intervention contamination instead of re-deriving frozen "
    "group identity from free-text tokens.\n"
    "2. Duration correction: publication states 21 weeks DIO-MASH diet before "
    "treatment plus another 6 weeks of treatment with diet feeding continued. "
    "Final HFD exposure under the prespecified start-to-collection definition is "
    "27 weeks. The 28-week value written in Batch06 v1 was an arithmetic/source "
    "interpretation error and is superseded by 27 weeks.\n"
    "\n"
    "No frozen biological eligibility, technical eligibility, sample membership, "
    "or contrast membership changed.\n"
)

# Status summary
status_rows = []
for gse, z in x.groupby("GSE"):
    sts = set(s for s in z["overall_adjudication_status"] if s)
    if any("PARTIAL_REVIEWED" in s for s in sts):
        st = "PARTIAL_REVIEWED"
    elif sts and all("ADJUDICATED" in s for s in sts):
        st = "FULLY_ADJUDICATED"
    else:
        st = "NOT_YET_REVIEWED"
    status_rows.append((gse, st, z["Contrast_id"].nunique()))

sdf = pd.DataFrame(status_rows, columns=["GSE", "status", "n_contrasts"])

print("===== HFD MODERATOR MANUAL ADJUDICATION BATCH06 v1.1 PASS =====")
print("Reviewed GSE this batch: 5")
print("Fully adjudicated this batch: 4")
print("Partially reviewed this batch: 1")
print("Cumulative fully adjudicated GSE:",
      int((sdf["status"] == "FULLY_ADJUDICATED").sum()))
print("Cumulative partial-reviewed GSE:",
      int((sdf["status"] == "PARTIAL_REVIEWED").sum()))
print("Cumulative not-yet-reviewed GSE:",
      int((sdf["status"] == "NOT_YET_REVIEWED").sum()))

print("\n===== BATCH06 FINAL STATUS =====")
cols = [
    "GSE", "Contrast_id", "overall_adjudication_status",
    "sex_final", "strain_final", "hfd_duration_weeks_final",
    "hfd_duration_weeks_final_confidence",
    "hfd_fat_kcal_percent_final", "control_diet_type_final",
]
print(b[cols].sort_values(["GSE", "Contrast_id"]).to_string(index=False))

print("\n===== SOURCE-INSUFFICIENT FIELDS =====")
print(pd.DataFrame(attempts).to_string(index=False))

print("\nKEY CORRECTIONS:")
print("- GSE268200: 21 weeks was pre-treatment HFD duration; publication states")
print("  another 6 weeks with diet continued, so final exposure = 27 weeks.")
print("  The earlier Batch06 v1 value of 28 weeks is superseded.")
print("- GSE285614: auto 8 weeks was age at HFD start; bulk liver RNA-seq was after")
print(f"  3 months HFD -> {duration_3mo} weeks (MEDIUM confidence).")
print("- GSE205846: TD.08806 is a defined 10.5 kcal% fat control, categorized as LFD.")
print("- GSE191004: SAFE 260 HF = 60% energy from fat; strain remains unresolved.")

print("\nEXPECTED DURATION EFFECT:")
print("- Strict-ready upgrades expected from GSE293219, GSE268200, GSE191004, GSE205846.")
print("- GSE285614 remains non-strict for duration because month->week conversion is MEDIUM.")
print("- Therefore strict-duration independent GSE should rise from 24 to 28 after recomputation.")

print("\nSaved:")
print(OUT.relative_to(ROOT))
print(ATTEMPT_OUT.relative_to(ROOT))
print(CORRECTION_NOTE.relative_to(ROOT))

print("\nIMPORTANT:")
print("- No strain/substrain was inferred where source mapping was insufficient.")
print("- No intervention/genotype arm was mixed into frozen primary contrasts.")
print("- No frozen biological/technical eligibility or contrast membership changed.")
