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

MASTER = (
    ROOT /
    "00_metadata/contrasts/moderators/manual_adjudication/"
    "HFD_primary105_moderator_manual_adjudication_batch03_v1.tsv"
)
GSM_AUTO = (
    ROOT /
    "00_metadata/contrasts/moderators/"
    "HFD_primary837_GSM_moderator_extraction_GEOauto_v1_1_1.tsv"
)

OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"
OUT = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch04_v1.tsv"
ATTEMPT_OUT = OUTDIR / "HFD_batch04_source_insufficient_review_log_v1.tsv"
DURATION_POLICY_OUT = OUTDIR / "HFD_duration_unit_standardization_policy_v1.txt"

if not MASTER.exists():
    raise FileNotFoundError(MASTER)
if not GSM_AUTO.exists():
    raise FileNotFoundError(GSM_AUTO)

x = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
g = pd.read_csv(GSM_AUTO, sep="\t", dtype=str).fillna("")

assert x["GSE"].nunique() == 64
assert x["Contrast_id"].nunique() == 105
assert g["GSM"].nunique() == 837
assert g["Contrast_id"].nunique() == 105

def fill_var(
    idx,
    var,
    manual_value,
    source_type,
    citation,
    evidence,
    confidence,
    status="MANUAL_CONFIRMED",
    final_value=None,
    final_note="",
):
    if final_value is None:
        final_value = manual_value

    x.loc[idx, f"{var}_manual_value"] = str(manual_value)
    x.loc[idx, f"{var}_manual_source_type"] = source_type
    x.loc[idx, f"{var}_manual_source_citation"] = citation
    x.loc[idx, f"{var}_manual_evidence"] = evidence
    x.loc[idx, f"{var}_manual_confidence"] = confidence
    x.loc[idx, f"{var}_adjudication_status"] = status
    x.loc[idx, f"{var}_final"] = str(final_value)
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
# Duration-unit policy
# ============================================================
DURATION_POLICY_OUT.write_text(
    "HFD duration unit standardization policy v1\n"
    "\n"
    "Definition:\n"
    "HFD duration is time from initiation of HFD feeding to tissue collection\n"
    "for the specific frozen HFD-vs-control contrast.\n"
    "\n"
    "Rules:\n"
    "1. Explicit weeks are retained directly.\n"
    "2. Explicit days are converted by days / 7.\n"
    "3. Calendar-month-only reporting is converted for quantitative meta-regression\n"
    "   using 365.25 / 12 / 7 = 4.3482142857 weeks per month.\n"
    "4. Month-derived values are marked MEDIUM confidence and the original\n"
    "   source wording is retained in manual evidence/final note.\n"
    "5. No conversion is used to infer an unreported duration.\n"
    "6. Sensitivity analysis may round month-derived durations to the nearest\n"
    "   whole week or exclude approximate-unit conversions.\n"
)

# ============================================================
# GSE145840 — COMPLETE
# Rajagopalan et al. J Clin Invest. 2020.
# PMID 32780721; DOI 10.1172/JCI137315
#
# Primary paper:
# - C57BL/6J male
# - HFD positive-control cohort from age 3 weeks
# - HFD feeding 14 weeks
# - chow comparator
#
# HFD fat composition is not stated in the accessible primary article text.
# A later curated study table explicitly identifies GSE145840 as
# 60 kcal% fat for 14 weeks; retain MEDIUM confidence for fat%.
# ============================================================
gse = "GSE145840"
idx = set_pub(
    gse,
    "32780721",
    "10.1172/JCI137315",
    "Metabolic effects of air pollution exposure and reversibility",
)
citation = (
    "Rajagopalan et al. J Clin Invest. 2020;130:6034-6040. "
    "PMID:32780721; DOI:10.1172/JCI137315"
)

t = selected_text(gse)
if not re.search(r"HFD|high[- ]fat", t, re.I):
    raise RuntimeError("GSE145840: HFD identity not detected in selected GSM metadata.")
if re.search(r"PM2\.?5|particulate|cessation", t, re.I):
    raise RuntimeError(
        "GSE145840: PM2.5/cessation samples detected in frozen primary contrast."
    )

fill_var(
    idx, "sex", "MALE",
    "PUBLICATION_METHODS", citation,
    "Primary publication identifies the early-life experiment as C57BL/6J male mice.",
    "HIGH",
)
fill_var(
    idx, "strain", "C57BL/6J",
    "PUBLICATION_METHODS", citation,
    "Primary publication explicitly reports C57BL/6J mice.",
    "HIGH",
)
fill_var(
    idx, "hfd_duration_weeks", "14",
    "PUBLICATION_METHODS", citation,
    "Three-week-old C57BL/6J mice were fed HFD for 14 weeks as the positive-control cohort.",
    "HIGH",
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60",
    "OTHER",
    (
        citation +
        "; Xu et al. Front Nutr. 2024;11:1381779. "
        "DOI:10.3389/fnut.2024.1381779, Table 1"
    ),
    "A later curated transcriptomic study table explicitly identifies GSE145840 "
    "as a 60 kcal% fat HFD model for 14 weeks. The accessible primary article "
    "confirms the 14-week HFD cohort but does not state the fat percentage.",
    "MEDIUM",
    final_note=(
        "Fat percentage supported by secondary curated dataset table rather than "
        "the accessible primary Methods text; retain for analysis with confidence flag."
    ),
)
fill_var(
    idx, "control_diet_type", "CHOW",
    "PUBLICATION_METHODS", citation,
    "Primary publication compares chow diet with HFD.",
    "HIGH",
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected GSMs contain the HFD-vs-diet-control stratum; PM2.5 and "
    "cessation arms are rejected by QA."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH04_ADJUDICATED"

# ============================================================
# GSE255354 — COMPLETE
# Zhu et al. Research Square preprint.
# DOI 10.21203/rs.3.rs-4130669/v1
#
# - 5-week-old C57BL/6 male mice
# - ND vs HFD, HFD = 60% calories from fat, XTHF60
# - first 9 weeks HFD, then vehicle-control HFD continues for 10 weeks
# - frozen primary HFD-vs-ND baseline therefore has 19 weeks total exposure
# - MEC1-treated HFD must not enter primary contrast
# ============================================================
gse = "GSE255354"
idx = set_pub(
    gse,
    "",
    "10.21203/rs.3.rs-4130669/v1",
    "Enhancement of ubiquitination-dependent mitophagy by unconventionally-acting PROTACs",
)
citation = (
    "Zhu et al. Research Square preprint. "
    "DOI:10.21203/rs.3.rs-4130669/v1; GEO:GSE255354"
)

t = selected_text(gse)
if not re.search(r"HFD|high[- ]fat", t, re.I):
    raise RuntimeError("GSE255354: HFD identity not detected in selected GSM metadata.")
if not re.search(r"\bND\b|normal diet|NCD", t, re.I):
    raise RuntimeError("GSE255354: normal-diet control not detected.")
if re.search(r"MEC1|MEC[-_ ]?1", t, re.I):
    raise RuntimeError(
        "GSE255354: MEC1-treated samples detected in frozen primary contrast."
    )

fill_var(
    idx, "sex", "MALE",
    "PUBLICATION_METHODS", citation,
    "DIO Methods specify 5-week-old C57BL/6 male mice.",
    "HIGH",
)
fill_var(
    idx, "strain", "C57BL/6",
    "PUBLICATION_METHODS", citation,
    "DIO Methods explicitly report C57BL/6; no J/N substrain is inferred.",
    "HIGH",
)
fill_var(
    idx, "hfd_duration_weeks", "19",
    "PUBLICATION_METHODS", citation,
    "Mice received ND or HFD for 9 weeks; HFD control mice then continued HFD "
    "during the subsequent 10-week vehicle/MEC1 phase. Frozen primary uses "
    "untreated/vehicle HFD control versus ND, giving 19 weeks HFD exposure.",
    "HIGH",
    final_note="9 weeks initial HFD + 10 weeks continued HFD = 19 weeks to tissue collection.",
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60",
    "PUBLICATION_METHODS", citation,
    "HFD was XTHF60 (Jiangshu Xietong), with 60% of calories from fat.",
    "HIGH",
)
fill_var(
    idx, "control_diet_type", "NORMAL_DIET",
    "PUBLICATION_METHODS", citation,
    "Comparator was explicitly termed normal diet (ND); Methods report 4.5% fat "
    "without a kcal unit, so that 4.5% value is not converted to kcal%.",
    "HIGH",
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Selected frozen contrast is ND vs HFD control; MEC1-treated HFD samples "
    "are explicitly rejected by QA."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH04_ADJUDICATED"

# ============================================================
# GSE285380 — COMPLETE with documented month->week conversion
# J Clin Invest. 2025;135:e193976.
# DOI 10.1172/JCI193976
#
# GEO GSM treatment protocol:
# - female C57BL/6
# - normal chow vs HFHSD
# - HFHSD = 29% Fat Kcal, Teklad
# - 3 months
# - frozen primary should be lean vs obese WITHOUT estrogen pellets
# ============================================================
gse = "GSE285380"
idx = set_pub(
    gse,
    "",
    "10.1172/JCI193976",
    "Estrogen and obesity synergistically suppress protein S via HIF1α, enhancing thrombosis potential",
)
citation = (
    "J Clin Invest. 2025;135:e193976. DOI:10.1172/JCI193976; "
    "GEO GSE285380 selected mouse liver GSM treatment protocol"
)

t = selected_text(gse)
if not re.search(r"Lean", t, re.I):
    raise RuntimeError("GSE285380: lean control identity not detected.")
if not re.search(r"Obese", t, re.I):
    raise RuntimeError("GSE285380: obese/HFD identity not detected.")
if re.search(r"Estrogen|Estradiol|\+P\b", t, re.I):
    raise RuntimeError(
        "GSE285380: estrogen-pelleted samples detected in frozen primary contrast."
    )

MONTH_TO_WEEK = 365.25 / 12.0 / 7.0
duration_3mo = 3.0 * MONTH_TO_WEEK
duration_3mo_s = f"{duration_3mo:.4f}"

fill_var(
    idx, "sex", "FEMALE",
    "GEO_GSM", citation,
    "Selected mouse liver GSM treatment protocol explicitly reports C57BL/6 female mice.",
    "HIGH",
)
fill_var(
    idx, "strain", "C57BL/6",
    "GEO_GSM", citation,
    "Selected mouse liver GSM treatment protocol explicitly reports C57BL/6; "
    "no J/N substrain is inferred.",
    "HIGH",
)
fill_var(
    idx, "hfd_duration_weeks", duration_3mo_s,
    "GEO_GSM", citation,
    "GEO reports normal chow or HFHSD feeding for 3 months. For quantitative "
    "meta-regression, 3 calendar months were standardized using "
    "365.25/12/7 weeks per month.",
    "MEDIUM",
    final_note=(
        f"Source unit='3 months'; standardized={duration_3mo_s} weeks. "
        "See HFD_duration_unit_standardization_policy_v1.txt. "
        "Sensitivity may round to 13 weeks or exclude approximate conversions."
    ),
)
fill_var(
    idx, "hfd_fat_kcal_percent", "29",
    "GEO_GSM", citation,
    "GEO treatment protocol explicitly reports HFHSD: 29% Fat Kcal, Teklad.",
    "HIGH",
)
fill_var(
    idx, "control_diet_type", "CHOW",
    "GEO_GSM", citation,
    "Selected lean comparator was explicitly fed normal chow diet.",
    "HIGH",
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected liver samples are lean/no-estrogen versus obese/no-estrogen; "
    "estrogen-pelleted cohorts are rejected by QA."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH04_ADJUDICATED"

# ============================================================
# GSE279693 — COMPLETE
# Liu et al. Metabolism. 2025;168:156258.
# PMID 40185397; DOI 10.1016/j.metabol.2025.156258
#
# - male C57BL/6J
# - D12492 60% fat
# - 16 weeks from age 8 weeks
# - standard normal chow
# - frozen primary must use sham-only NS vs HS, excluding IRI NI/HI
# ============================================================
gse = "GSE279693"
idx = set_pub(
    gse,
    "40185397",
    "10.1016/j.metabol.2025.156258",
    "Hexosamine biosynthesis dysfunction-induced LIFR N-glycosylation deficiency exacerbates steatotic liver ischemia/reperfusion injury",
)
citation = (
    "Liu et al. Metabolism. 2025;168:156258. "
    "PMID:40185397; DOI:10.1016/j.metabol.2025.156258"
)

t = selected_text(gse)
# Expected clean sham baseline sample titles are NS1-3 and HS1-3.
if not re.search(r"\bNS[123]\b", t, re.I):
    raise RuntimeError("GSE279693: normal-diet sham (NS) samples not detected.")
if not re.search(r"\bHS[123]\b", t, re.I):
    raise RuntimeError("GSE279693: HFD sham (HS) samples not detected.")
if re.search(r"\bNI[123]\b|\bHI[123]\b", t, re.I):
    raise RuntimeError(
        "GSE279693: ischemia-reperfusion (NI/HI) samples detected in frozen primary contrast."
    )

fill_var(
    idx, "sex", "MALE",
    "PUBLICATION_METHODS", citation,
    "MASLD model Methods explicitly report male mice.",
    "HIGH",
)
fill_var(
    idx, "strain", "C57BL/6J",
    "PUBLICATION_METHODS", citation,
    "All experimental mouse lines are reported on a C57BL/6J background.",
    "HIGH",
)
fill_var(
    idx, "hfd_duration_weeks", "16",
    "PUBLICATION_METHODS", citation,
    "Male mice were fed HFD continuously for 16 weeks starting at 8 weeks of age.",
    "HIGH",
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60",
    "PUBLICATION_METHODS", citation,
    "HFD was Research Diets D12492, reported as 60% fat.",
    "HIGH",
)
fill_var(
    idx, "control_diet_type", "CHOW",
    "PUBLICATION_METHODS", citation,
    "Control animals received a standard normal chow diet.",
    "HIGH",
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen primary contrast verified NS (normal-diet sham) vs HS (HFD sham); "
    "NI/HI ischemia-reperfusion samples are explicitly rejected by QA."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH04_ADJUDICATED"

# ============================================================
# GSE132056 — PARTIAL
# Chaurasia et al. Science. 2019.
# PMID 31273070; DOI 10.1126/science.aav3722
#
# Frozen primary: NCD control vs HFD control, not liver-specific Degs1 KO.
# - mice backcrossed 10x onto C57BL/6 background
# - HFD from 4 weeks of age; euthanasia/tissue harvest at 16 weeks
#   => 12 weeks HFD exposure
# - sex and HFD fat composition not explicitly supported in sources reviewed
# ============================================================
gse = "GSE132056"
idx = set_pub(
    gse,
    "31273070",
    "10.1126/science.aav3722",
    "Targeting a ceramide double bond improves insulin resistance and hepatic steatosis",
)
citation = (
    "Chaurasia et al. Science. 2019;365:386-392. "
    "PMID:31273070; DOI:10.1126/science.aav3722; GEO:GSE132056"
)

t = selected_text(gse)
if not re.search(r"HFD control", t, re.I):
    raise RuntimeError("GSE132056: HFD control samples not detected.")
if not re.search(r"NCD control", t, re.I):
    raise RuntimeError("GSE132056: NCD control samples not detected.")
if re.search(r"knockout", t, re.I):
    raise RuntimeError(
        "GSE132056: Degs1 knockout sample detected in frozen primary contrast."
    )

fill_var(
    idx, "strain", "C57BL/6",
    "PUBLICATION_METHODS", citation,
    "Conditional Degs1 mice were backcrossed 10 generations onto the C57BL/6 background; "
    "no J/N substrain is inferred.",
    "HIGH",
)
fill_var(
    idx, "hfd_duration_weeks", "12",
    "PUBLICATION_METHODS", citation,
    "The tissue-specific HFD cohort was placed on HFD at 4 weeks of age and "
    "euthanized/tissues harvested at 16 weeks of age.",
    "HIGH",
)
fill_var(
    idx, "control_diet_type", "CHOW",
    "GEO_GSM", "GSE132056 selected GSM characteristics + primary publication",
    "Frozen comparator GSMs are explicitly labeled normal chow diet (NCD) control.",
    "HIGH",
)
mark_not_reported(
    idx, "sex",
    "PUBLICATION_METHODS", citation,
    "Sex of the specific frozen RNA-seq NCD-control/HFD-control cohort was not "
    "explicitly supported by the sources reviewed; no inference is made."
)
mark_not_reported(
    idx, "hfd_fat_kcal_percent",
    "PUBLICATION_METHODS", citation,
    "The reviewed primary article/GEO records identify an obesogenic HFD but do "
    "not explicitly provide kcal% fat for this frozen RNA-seq cohort."
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected samples verified as NCD control vs HFD control; "
    "Degs1 knockout samples are excluded by QA."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH04_PARTIAL_REVIEWED"

attempts.append({
    "GSE": "GSE132056",
    "unresolved_fields": "sex,hfd_fat_kcal_percent",
    "reason": (
        "Specific frozen RNA-seq cohort sex and kcal% fat were not explicitly "
        "supported by the reviewed primary/GEO sources; no inference made."
    ),
})

# ------------------------------------------------------------
# QA
# ------------------------------------------------------------
assert x["GSE"].nunique() == 64
assert x["Contrast_id"].nunique() == 105

batch_gses = [
    "GSE145840",
    "GSE255354",
    "GSE285380",
    "GSE279693",
    "GSE132056",
]
b = x[x["GSE"].isin(batch_gses)].copy()

assert b["GSE"].nunique() == 5
assert b["Contrast_id"].nunique() == 5

# Four complete cases.
for gse in ["GSE145840", "GSE255354", "GSE285380", "GSE279693"]:
    z = b[b["GSE"].eq(gse)]
    for var in [
        "sex", "strain", "hfd_duration_weeks",
        "hfd_fat_kcal_percent", "control_diet_type"
    ]:
        assert (z[f"{var}_final"] != "").all(), (gse, var)
    assert set(z["overall_adjudication_status"]) == {"BATCH04_ADJUDICATED"}

# Partial case.
z = b[b["GSE"].eq("GSE132056")]
assert set(z["strain_final"]) == {"C57BL/6"}
assert set(z["hfd_duration_weeks_final"]) == {"12"}
assert set(z["control_diet_type_final"]) == {"CHOW"}
assert set(z["sex_final"]) == {""}
assert set(z["hfd_fat_kcal_percent_final"]) == {""}
assert set(z["overall_adjudication_status"]) == {"BATCH04_PARTIAL_REVIEWED"}

# Specific values.
z = b[b["GSE"].eq("GSE145840")]
assert set(z["sex_final"]) == {"MALE"}
assert set(z["strain_final"]) == {"C57BL/6J"}
assert set(z["hfd_duration_weeks_final"]) == {"14"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"60"}

z = b[b["GSE"].eq("GSE255354")]
assert set(z["sex_final"]) == {"MALE"}
assert set(z["strain_final"]) == {"C57BL/6"}
assert set(z["hfd_duration_weeks_final"]) == {"19"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"60"}

z = b[b["GSE"].eq("GSE285380")]
assert set(z["sex_final"]) == {"FEMALE"}
assert set(z["strain_final"]) == {"C57BL/6"}
assert set(z["hfd_duration_weeks_final"]) == {duration_3mo_s}
assert set(z["hfd_fat_kcal_percent_final"]) == {"29"}

z = b[b["GSE"].eq("GSE279693")]
assert set(z["sex_final"]) == {"MALE"}
assert set(z["strain_final"]) == {"C57BL/6J"}
assert set(z["hfd_duration_weeks_final"]) == {"16"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"60"}
assert set(z["control_diet_type_final"]) == {"CHOW"}

# Preserve previously-reviewed GSE statuses.
prev = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
for gse in [
    "GSE108688", "GSE236201", "GSE278920",
    "GSE218932", "GSE246221", "GSE294746", "GSE296500",
    "GSE198686", "GSE205021", "GSE267072", "GSE284023", "GSE304736",
]:
    a = prev[prev["GSE"].eq(gse)].sort_values("Contrast_id")
    bb = x[x["GSE"].eq(gse)].sort_values("Contrast_id")
    assert a["overall_adjudication_status"].tolist() == bb["overall_adjudication_status"].tolist(), gse

x.to_csv(OUT, sep="\t", index=False)
pd.DataFrame(attempts).to_csv(ATTEMPT_OUT, sep="\t", index=False)

# Current status counts after Batch04.
status_by_gse = []
for gse, z in x.groupby("GSE"):
    sts = set(s for s in z["overall_adjudication_status"] if s)
    if any("PARTIAL_REVIEWED" in s for s in sts):
        st = "PARTIAL_REVIEWED"
    elif sts and all("ADJUDICATED" in s for s in sts):
        st = "FULLY_ADJUDICATED"
    else:
        st = "NOT_YET_REVIEWED"
    status_by_gse.append((gse, st, z["Contrast_id"].nunique()))

sdf = pd.DataFrame(status_by_gse, columns=["GSE","status","n_contrasts"])

print("===== HFD MODERATOR MANUAL ADJUDICATION BATCH04 PASS =====")
print("Reviewed GSE this batch: 5")
print("Fully adjudicated this batch: 4")
print("Partially reviewed this batch: 1")
print("Cumulative fully adjudicated GSE:",
      int((sdf["status"] == "FULLY_ADJUDICATED").sum()))
print("Cumulative partial-reviewed GSE:",
      int((sdf["status"] == "PARTIAL_REVIEWED").sum()))
print("Cumulative not-yet-reviewed GSE:",
      int((sdf["status"] == "NOT_YET_REVIEWED").sum()))

print("\n===== BATCH04 FINAL STATUS =====")
cols = [
    "GSE", "Contrast_id", "overall_adjudication_status",
    "sex_final", "strain_final", "hfd_duration_weeks_final",
    "hfd_fat_kcal_percent_final", "control_diet_type_final",
]
print(b[cols].sort_values(["GSE","Contrast_id"]).to_string(index=False))

print("\n===== SOURCE-INSUFFICIENT FIELDS =====")
print(pd.DataFrame(attempts).to_string(index=False))

print("\n===== DURATION UNIT POLICY =====")
print(f"weeks_per_calendar_month = {MONTH_TO_WEEK:.10f}")
print(f"GSE285380 3 months -> {duration_3mo_s} weeks (MEDIUM confidence)")

print("\nSaved:")
print(OUT.relative_to(ROOT))
print(ATTEMPT_OUT.relative_to(ROOT))
print(DURATION_POLICY_OUT.relative_to(ROOT))

print("\nIMPORTANT:")
print("- GSE145840 fat%=60 has MEDIUM confidence because it is supported by")
print("  a later curated dataset table; primary paper directly confirms sex, strain,")
print("  and 14-week HFD but not the fat percentage in accessible Methods text.")
print("- GSE285380 month-to-week conversion is explicit and sensitivity-testable.")
print("- GSE132056 sex/fat remain unresolved.")
print("- No J/N substrain was inferred when source reports only C57BL/6.")
print("- No frozen eligibility or contrast membership changed.")
