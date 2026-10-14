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
    "HFD_primary105_moderator_manual_adjudication_batch02_v1.tsv"
)
GSM_AUTO = (
    ROOT /
    "00_metadata/contrasts/moderators/"
    "HFD_primary837_GSM_moderator_extraction_GEOauto_v1_1_1.tsv"
)
OUT = (
    ROOT /
    "00_metadata/contrasts/moderators/manual_adjudication/"
    "HFD_primary105_moderator_manual_adjudication_batch03_v1.tsv"
)
ATTEMPT_OUT = (
    ROOT /
    "00_metadata/contrasts/moderators/manual_adjudication/"
    "HFD_batch03_source_insufficient_review_log_v1.tsv"
)

if not MASTER.exists():
    raise FileNotFoundError(MASTER)
if not GSM_AUTO.exists():
    raise FileNotFoundError(GSM_AUTO)

x = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
g = pd.read_csv(GSM_AUTO, sep="\t", dtype=str).fillna("")

assert x["Contrast_id"].nunique() == 105
assert x["GSE"].nunique() == 64
assert g["GSM"].nunique() == 837

def fill_var(
    idx,
    var,
    manual_value,
    source_type,
    citation,
    evidence,
    confidence,
    status,
    final_note="",
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
    # Intentionally keep *_final blank so future coverage logic does not
    # count this moderator as analysis-ready.
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

def selected_text(gse, cid=None):
    z = g[g["GSE"].eq(gse)].copy()
    if cid is not None:
        z = z[z["Contrast_id"].eq(cid)]
    vals = []
    for c in ["sample_title_raw", "sample_characteristics_raw"]:
        if c in z.columns:
            vals.extend(v for v in z[c].tolist() if v)
    return " || ".join(vals)

attempts = []

# ============================================================
# GSE198686 — COMPLETE
# Sci Adv 2023; PMID 37672589
# RNA-seq cohort explicitly all male.
# C57BL/6N Nme2 model.
# HFD SAFE U8978 v19 = SAFE 260HF = 60% energy from fat.
# ND SAFE A03 = standard chow/home-cage diet.
# HFD duration = 6 weeks.
# ============================================================
gse = "GSE198686"
idx = set_pub(
    gse,
    "37672589",
    "10.1126/sciadv.adh0140",
    "Nucleoside diphosphate kinases 1 and 2 regulate a protective liver response to a high-fat diet",
)
citation = (
    "Iuso et al. Sci Adv. 2023;9:eadh0140. "
    "PMID:37672589; DOI:10.1126/sciadv.adh0140"
)

# Verify frozen selected contrast is WT/normal vs WT/HFD, not KO.
t = selected_text(gse)
if re.search(r"nme2ko|Nme2[-_/ ]?KO|knockout", t, re.I):
    raise RuntimeError(
        "GSE198686 frozen primary contrast appears to contain Nme2 KO samples; "
        "manual review required before applying WT-only adjudication."
    )
if not re.search(r"nme2wt|wild.?type|\bWT\b", t, re.I):
    raise RuntimeError(
        "GSE198686 WT identity not detectable in selected GSM metadata."
    )

fill_var(
    idx, "sex", "MALE",
    "PUBLICATION_METHODS", citation,
    "Quantification/statistical section states RNA-seq used four male mice "
    "per genotype/diet condition.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "strain", "C57BL/6N",
    "PUBLICATION_METHODS", citation,
    "Nme2 model: C57BL/6N-Nme2tm1.1(KOMP)Vlcg/MbpMmucd; "
    "frozen primary selected samples are WT.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "hfd_duration_weeks", "6",
    "PUBLICATION_METHODS", citation,
    "WT and Nme2 KO mice were fed ND or HFD for 6 weeks; RNA-seq was "
    "performed after the 6-week challenge.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60",
    "DIET_MANUFACTURER_CODE",
    citation + "; SAFE U8978 Version 19 / SAFE 260HF",
    "Publication specifies SAFE U8978 version 19; manufacturer identifies "
    "SAFE 260HF U8978 v19 as 60% of energy from fat.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "control_diet_type", "CHOW",
    "PUBLICATION_METHODS", citation,
    "Normal diet was SAFE A03, a standard rodent breeding/home-cage diet.",
    "HIGH", "MANUAL_CONFIRMED",
)
x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected GSM metadata verified WT normal-diet vs WT HFD; "
    "Nme2 KO strata are not used in this primary contrast."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH03_ADJUDICATED"

# ============================================================
# GSE284023 — COMPLETE
# Cell Metab 2026; PMID 41421349
# male WT C57BL/6J
# STD TD2914 = 13% energy fat
# HFD TD.06414 = 60.3% energy fat
# short-term diet = 1 week; long-term = 18 weeks
# Frozen selected contrasts must be untreated STD-vs-HFD only.
# ============================================================
gse = "GSE284023"
idx = set_pub(
    gse,
    "41421349",
    "10.1016/j.cmet.2025.11.012",
    "Enhanced non-enzymatic H2S generation extends lifespan and healthspan in male mice",
)
citation = (
    "Cáliz-Molina et al. Cell Metab. 2026;38:350-369.e8. "
    "PMID:41421349; DOI:10.1016/j.cmet.2025.11.012"
)

fill_var(
    idx, "sex", "MALE",
    "GEO_GSM", "GSE284023 growth protocol + publication",
    "GEO growth protocol identifies the RNA-seq diet cohorts as male mice; "
    "publication title/Methods are concordant.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "strain", "C57BL/6J",
    "PUBLICATION_METHODS", citation,
    "Wild-type C57BL/6J mice were purchased from Charles River.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60.3",
    "PUBLICATION_METHODS", citation,
    "HFD TD.06414 contained 60.3% of calories from fat.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "control_diet_type", "CHOW",
    "PUBLICATION_METHODS", citation,
    "STD was Envigo TD2914 standard rodent chow, 13% calories from fat.",
    "HIGH", "MANUAL_CONFIRMED",
)

cids = sorted(x.loc[idx, "Contrast_id"].unique())
assert len(cids) == 2, cids

for cid in cids:
    t = selected_text(gse, cid)

    # Guard against treatment/sugar arms.
    if re.search(r"\bDAS\b|\bDAD\b|\bDAT\b|\bHFDHS\b|\bHS\b", t, re.I):
        raise RuntimeError(
            f"{cid}: intervention or high-sucrose sample detected in frozen "
            "primary contrast; manual inspection required."
        )
    if not re.search(r"\bHFD\b", t, re.I):
        raise RuntimeError(f"{cid}: HFD sample not detectable.")
    if not re.search(r"\bSTD\b", t, re.I):
        raise RuntimeError(f"{cid}: STD control sample not detectable.")

    has_short = bool(re.search(r"short.?term", t, re.I))
    has_long = bool(re.search(r"long.?term", t, re.I))

    if has_short and not has_long:
        duration = "1"
        note = "Selected frozen GSMs are untreated STD-vs-HFD short-term samples."
    elif has_long and not has_short:
        duration = "18"
        note = "Selected frozen GSMs are untreated STD-vs-HFD long-term samples."
    else:
        raise RuntimeError(
            f"{cid}: cannot uniquely resolve short-term vs long-term selected GSMs."
        )

    fill_var(
        x["Contrast_id"].eq(cid),
        "hfd_duration_weeks",
        duration,
        "GEO_GSM",
        "GSE284023 selected GSM growth protocol + publication",
        note + " GEO growth protocol states 1 week for short-term and "
        "18 weeks for long-term diet cohorts.",
        "HIGH",
        "PAPER_GEO_CONFLICT_RESOLVED",
        final_note=(
            "Prior GEO-auto conflict mixed age/treatment text. Duration is "
            "resolved from the selected untreated STD/HFD cohort label."
        ),
    )

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Selected frozen contrasts are untreated STD-vs-HFD only; DAD/DAT/DAS "
    "and high-sucrose/HFDHS groups are guarded against by script."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH03_ADJUDICATED"

# ============================================================
# GSE304736 — PARTIAL (strain not explicitly reported)
# Metabolism 2025; PMID 40446956
#
# male, HFD D12492 60 kcal% fat, standard chow NIH-31.
# For the obese RNA-seq HFD-vs-STD vehicle comparison, tissue collection
# occurs after 14 weeks HFD + 10 days continued HFD during vehicle phase:
# 108 days = 15.428571 weeks under our prespecified duration definition.
# ============================================================
gse = "GSE304736"
idx = set_pub(
    gse,
    "40446956",
    "10.1016/j.metabol.2025.156308",
    "Hepatic CB1 receptor signaling triggers Gi/oα-mediated lipolysis in lean mice but Gsα-mediated lipogenesis in obese mice",
)
citation = (
    "Liu et al. Metabolism. 2025;170:156308. "
    "PMID:40446956; DOI:10.1016/j.metabol.2025.156308"
)

fill_var(
    idx, "sex", "MALE",
    "PUBLICATION_METHODS", citation,
    "Male six-week-old mice were placed on standard chow or HFD.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "hfd_duration_weeks", "15.428571",
    "PUBLICATION_METHODS", citation,
    "Mice received HFD for 14 weeks, then vehicle/JD5037 for an additional "
    "10 days with continued HFD. Frozen HFD-vs-STD RNA-seq comparison is "
    "required to be the vehicle baseline; total HFD exposure to collection "
    "= 14 weeks + 10 days = 108 days = 15.428571 weeks.",
    "HIGH", "MANUAL_CONFIRMED",
    final_note=(
        "Duration follows prespecified definition: HFD initiation to tissue "
        "collection, not merely pre-treatment HFD phase."
    ),
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60",
    "PUBLICATION_METHODS", citation,
    "HFD was Research Diets D12492, 60% kcal from fat.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "control_diet_type", "CHOW",
    "PUBLICATION_METHODS", citation,
    "STD was NIH-31 rodent standard chow.",
    "HIGH", "MANUAL_CONFIRMED",
)
mark_not_reported(
    idx, "strain",
    "PUBLICATION_METHODS", citation,
    "The publication describes hCB1R knockout/floxed littermates but does "
    "not explicitly report a C57BL/6 J/N substrain for the RNA-seq cohort. "
    "No substrain is inferred from transgenic components."
)

# Guard against drug-treated contrast.
t = selected_text(gse)
if re.search(r"JD5037|\bJD\b|MRI-1891|MRI1891", t, re.I):
    raise RuntimeError(
        "GSE304736 frozen primary contrast contains CB1R-antagonist-treated "
        "samples; expected untreated/vehicle HFD-vs-STD baseline."
    )

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Publication defines STD/HFD RNA-seq baseline; selected frozen metadata "
    "must not contain JD5037/MRI-1891 intervention samples."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH03_PARTIAL_REVIEWED"

# ============================================================
# GSE205021 — PARTIAL
# No linked primary publication in GEO.
# GEO explicitly supports 14-week HFD vs normal-diet control.
# Sex, strain, HFD fat composition are not reported in the public record
# we can support, so they remain unresolved.
# ============================================================
gse = "GSE205021"
idx = x["GSE"].eq(gse)

fill_var(
    idx, "hfd_duration_weeks", "14",
    "GEO_SERIES", "GSE205021 GEO Series overall design",
    "High-fat and fenofibrate groups received a 14-week HFD intervention; "
    "control group received normal diet for 14 weeks. Frozen primary uses "
    "HFD vs control, not HFD+fenofibrate.",
    "HIGH", "PAPER_GEO_CONFLICT_RESOLVED",
    final_note=(
        "Prior GEO-auto conflict captured the final 4-week fenofibrate phase; "
        "the untreated HFD-vs-control exposure is 14 weeks."
    ),
)
fill_var(
    idx, "control_diet_type", "NORMAL_DIET",
    "GEO_SERIES", "GSE205021 GEO Series overall design",
    "Control group is explicitly described as normal diet; composition/code "
    "is not reported.",
    "MEDIUM", "MANUAL_CONFIRMED",
)
for var in ["sex", "strain", "hfd_fat_kcal_percent"]:
    mark_not_reported(
        idx, var,
        "GEO_SERIES", "GSE205021 GEO Series/GSM public metadata",
        f"{var} was not explicitly reported in the public GEO Series/GSM "
        "metadata and no linked primary publication was available for "
        "source-supported adjudication."
    )

t = selected_text(gse)
if re.search(r"fenofibrate|\bFen\b", t, re.I):
    raise RuntimeError(
        "GSE205021 frozen primary contrast contains fenofibrate-treated sample."
    )

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected samples are Control vs High fat diet; fenofibrate arm "
    "is excluded from the primary contrast."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH03_PARTIAL_REVIEWED"

# ============================================================
# GSE267072 — REVIEWED, SOURCE INSUFFICIENT FOR MAIN MODERATORS
# GEO supports NCD-control vs HFD-control selected samples, but public GEO
# record does not report sex, strain, duration, or fat kcal%.
# Do not infer them from older Falcor papers or unrelated cohorts.
# ============================================================
gse = "GSE267072"
idx = x["GSE"].eq(gse)

# Verify selected primary uses Ctrl genotype, not Falcor KO.
t = selected_text(gse)
if re.search(r"Falcor.?KO|HFD_KO|NCD_KO", t, re.I):
    raise RuntimeError(
        "GSE267072 frozen primary contrast appears to include Falcor-KO samples; "
        "expected Ctrl NCD-vs-HFD baseline."
    )
if not re.search(r"HFD_Ctrl", t, re.I) or not re.search(r"NCD_Ctrl", t, re.I):
    raise RuntimeError(
        "GSE267072 selected Ctrl NCD/HFD identities not detectable."
    )

fill_var(
    idx, "control_diet_type", "NORMAL_DIET",
    "GEO_GSM", "GSE267072 selected GSM titles",
    "Selected control samples are explicitly labeled Liver_NCD_Ctrl; "
    "the public record does not provide a diet code/composition.",
    "MEDIUM", "MANUAL_CONFIRMED",
    final_note="NCD retained as normal-diet control; not forced to CHOW vs LFD."
)
for var in ["sex", "strain", "hfd_duration_weeks", "hfd_fat_kcal_percent"]:
    mark_not_reported(
        idx, var,
        "GEO_SERIES", "GSE267072 GEO Series/GSM public metadata",
        f"{var} is not explicitly reported in the public GEO record. "
        "No linked primary publication is listed in GEO; no inference made."
    )

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Selected frozen samples verified as NCD_Ctrl vs HFD_Ctrl; Falcor-KO "
    "samples are excluded from the primary contrast."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH03_PARTIAL_REVIEWED"

attempts.extend([
    {
        "GSE": "GSE304736",
        "unresolved_fields": "strain",
        "reason": "Publication does not explicitly report RNA-seq cohort substrain; no J/N inference."
    },
    {
        "GSE": "GSE205021",
        "unresolved_fields": "sex,strain,hfd_fat_kcal_percent",
        "reason": "Public GEO record lacks these fields; no linked primary publication."
    },
    {
        "GSE": "GSE267072",
        "unresolved_fields": "sex,strain,hfd_duration_weeks,hfd_fat_kcal_percent",
        "reason": "Public GEO record lacks these fields; GEO lists no primary publication."
    },
])

# ------------------------------------------------------------
# QA
# ------------------------------------------------------------
assert x["Contrast_id"].nunique() == 105
assert x["GSE"].nunique() == 64

# Complete cases
for gse in ["GSE198686", "GSE284023"]:
    z = x[x["GSE"].eq(gse)]
    for var in [
        "sex", "strain", "hfd_duration_weeks",
        "hfd_fat_kcal_percent", "control_diet_type"
    ]:
        assert (z[f"{var}_final"] != "").all(), (gse, var)

# Specific checks
z = x[x["GSE"].eq("GSE198686")]
assert set(z["sex_final"]) == {"MALE"}
assert set(z["strain_final"]) == {"C57BL/6N"}
assert set(z["hfd_duration_weeks_final"]) == {"6"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"60"}
assert set(z["control_diet_type_final"]) == {"CHOW"}

z = x[x["GSE"].eq("GSE284023")]
assert set(z["sex_final"]) == {"MALE"}
assert set(z["strain_final"]) == {"C57BL/6J"}
assert set(z["hfd_duration_weeks_final"]) == {"1", "18"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"60.3"}
assert set(z["control_diet_type_final"]) == {"CHOW"}

z = x[x["GSE"].eq("GSE304736")]
assert set(z["sex_final"]) == {"MALE"}
assert set(z["hfd_duration_weeks_final"]) == {"15.428571"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"60"}
assert set(z["control_diet_type_final"]) == {"CHOW"}
assert set(z["strain_final"]) == {""}

z = x[x["GSE"].eq("GSE205021")]
assert set(z["hfd_duration_weeks_final"]) == {"14"}
assert set(z["control_diet_type_final"]) == {"NORMAL_DIET"}
assert set(z["sex_final"]) == {""}
assert set(z["strain_final"]) == {""}
assert set(z["hfd_fat_kcal_percent_final"]) == {""}

x.to_csv(OUT, sep="\t", index=False)
pd.DataFrame(attempts).to_csv(ATTEMPT_OUT, sep="\t", index=False)

print("===== HFD MODERATOR MANUAL ADJUDICATION BATCH03 PASS =====")
print("Reviewed GSE: 5")
print("Fully adjudicated GSE this batch: 2")
print("Partially reviewed GSE this batch: 3")

print("\n===== BATCH03 STATUS =====")
cols = [
    "GSE", "Contrast_id", "overall_adjudication_status",
    "sex_final", "strain_final", "hfd_duration_weeks_final",
    "hfd_fat_kcal_percent_final", "control_diet_type_final"
]
print(
    x[x["GSE"].isin([
        "GSE198686","GSE284023","GSE304736","GSE205021","GSE267072"
    ])][cols].sort_values(["GSE","Contrast_id"]).to_string(index=False)
)

print("\n===== SOURCE-INSUFFICIENT FIELDS =====")
print(pd.DataFrame(attempts).to_string(index=False))

print("\nSaved:")
print(OUT.relative_to(ROOT))
print(ATTEMPT_OUT.relative_to(ROOT))

print("\nIMPORTANT:")
print("- Partial review is intentionally not counted as complete adjudication.")
print("- Blank final fields remain non-analysis-ready.")
print("- No missing sex/strain/diet composition was inferred.")
print("- No frozen biological/technical eligibility or contrast membership changed.")
