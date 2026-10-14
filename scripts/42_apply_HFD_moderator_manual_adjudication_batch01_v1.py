#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=__import__("os").environ.get("HFD_PROJECT_ROOT", ".")
)
args = ap.parse_args()

ROOT = Path(args.root)
INDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"

TEMPLATE = INDIR / "HFD_primary105_moderator_manual_adjudication_template_v1.tsv"
OUT = INDIR / "HFD_primary105_moderator_manual_adjudication_batch01_v1.tsv"
VOCAB = INDIR / "HFD_moderator_manual_adjudication_vocab_v1.tsv"
VOCAB_OUT = INDIR / "HFD_moderator_manual_adjudication_vocab_v1_1.tsv"

if not TEMPLATE.exists():
    raise FileNotFoundError(TEMPLATE)

x = pd.read_csv(TEMPLATE, sep="\t", dtype=str).fillna("")

assert x["GSE"].nunique() == 64
assert x["Contrast_id"].nunique() == 105

# ------------------------------------------------------------
# Controlled-vocabulary extension:
# Distinguish parser/auto-extraction error from true source conflict.
# ------------------------------------------------------------
if VOCAB.exists():
    v = pd.read_csv(VOCAB, sep="\t", dtype=str).fillna("")
    new = pd.DataFrame([{
        "category": "adjudication_status",
        "allowed_value": "AUTO_EXTRACTION_ERROR_CORRECTED"
    }])
    v = pd.concat([v, new], ignore_index=True).drop_duplicates()
    v.to_csv(VOCAB_OUT, sep="\t", index=False)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def fill_var(
    idx,
    var,
    manual_value,
    source_type,
    citation,
    evidence,
    manual_confidence,
    adjudication_status,
    final_value=None,
    final_confidence=None,
    final_source=None,
    final_note="",
):
    if final_value is None:
        final_value = manual_value
    if final_confidence is None:
        final_confidence = manual_confidence
    if final_source is None:
        final_source = source_type

    x.loc[idx, f"{var}_manual_value"] = str(manual_value)
    x.loc[idx, f"{var}_manual_source_type"] = source_type
    x.loc[idx, f"{var}_manual_source_citation"] = citation
    x.loc[idx, f"{var}_manual_evidence"] = evidence
    x.loc[idx, f"{var}_manual_confidence"] = manual_confidence
    x.loc[idx, f"{var}_adjudication_status"] = adjudication_status
    x.loc[idx, f"{var}_final"] = str(final_value)
    x.loc[idx, f"{var}_final_confidence"] = final_confidence
    x.loc[idx, f"{var}_final_source"] = final_source
    x.loc[idx, f"{var}_final_note"] = final_note

def set_pub(gse, pmid, doi, title):
    idx = x["GSE"].eq(gse)
    x.loc[idx, "publication_pmid"] = pmid
    x.loc[idx, "publication_doi"] = doi
    x.loc[idx, "publication_title"] = title
    return idx

# ============================================================
# BATCH 01
# GSE108688, GSE278920, GSE236201
# ============================================================

# ------------------------------------------------------------
# GSE108688
# Quagliarini et al., Molecular Cell 2019
# PMID 31706703; DOI 10.1016/j.molcel.2019.10.007
#
# Publication Methods:
# - male C57BL/6J
# - control D12329, 11 kcal% fat w/sucrose
# - HFD D12331, 58 kcal% fat w/sucrose
# - diet duration 12 weeks
#
# Important: GEO-auto HFD fat=11 was an extraction error caused by
# picking the control-diet fat percentage from shared/series text.
# This is NOT a GEO-vs-paper biological metadata conflict.
# ------------------------------------------------------------
gse = "GSE108688"
idx = set_pub(
    gse,
    "31706703",
    "10.1016/j.molcel.2019.10.007",
    "Cistromic Reprogramming of the Diurnal Glucocorticoid Hormone Response by High-Fat Diet",
)

citation = (
    "Quagliarini et al. Mol Cell. 2019;76:531-545.e5. "
    "PMID:31706703; DOI:10.1016/j.molcel.2019.10.007"
)

fill_var(
    idx, "sex", "MALE",
    "PUBLICATION_METHODS", citation,
    "Animal Experiments: male mice were used for the HFD/control liver study.",
    "HIGH", "MANUAL_CONFIRMED",
)

fill_var(
    idx, "hfd_duration_weeks", "12",
    "PUBLICATION_METHODS", citation,
    "Animals were fed control diet or HFD for 12 weeks.",
    "HIGH", "CONSISTENT",
)

fill_var(
    idx, "hfd_fat_kcal_percent", "58",
    "PUBLICATION_METHODS", citation,
    "HFD: 58 kcal% fat with sucrose, Research Diets D12331. "
    "Control diet D12329 contains 11 kcal% fat.",
    "HIGH", "AUTO_EXTRACTION_ERROR_CORRECTED",
    final_note=(
        "GEO-auto value 11% represented CONTROL diet D12329, not HFD. "
        "HFD D12331 is 58 kcal% fat."
    ),
)

fill_var(
    idx, "control_diet_type", "LFD",
    "PUBLICATION_METHODS", citation,
    "Control diet: Research Diets D12329, 11 kcal% fat with sucrose.",
    "HIGH", "MANUAL_CONFIRMED",
    final_note="Purified low-fat/control diet; not standard chow.",
)

# C01-C05 are WT/littermate-control strata from the Alb-Cre/GR-flox arm.
# GEO GSMs lack strain field, while the paper's animal section/model table
# supports C57BL/6J context. Keep confidence MEDIUM rather than silently
# treating it as direct GSM-level evidence.
idx_108_first = idx & x["Contrast_id"].isin(
    [f"GSE108688__C{i:02d}" for i in range(1, 6)]
)
fill_var(
    idx_108_first, "strain", "C57BL/6J",
    "PUBLICATION_METHODS", citation,
    "Publication Animal Experiments/Experimental Models report C57BL/6J; "
    "selected WT/littermate-control GSMs do not themselves carry a strain field.",
    "MEDIUM", "MANUAL_CONFIRMED",
    final_note="Publication-level confirmation; GSM-level strain absent.",
)

# C06-C10 have direct GSM-level C57BL/6J metadata.
idx_108_second = idx & x["Contrast_id"].isin(
    [f"GSE108688__C{i:02d}" for i in range(6, 11)]
)
fill_var(
    idx_108_second, "strain", "C57BL/6J",
    "GEO_GSM", "GSE108688 selected GSM characteristics",
    "Selected GSM characteristics explicitly report strain: C57BL/6J.",
    "HIGH", "CONSISTENT",
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Selected GSM titles/characteristics distinguish HFD D12331 from "
    "control/LFD D12329 and match publication design."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH01_ADJUDICATED"
x.loc[idx, "general_note"] = (
    "Circadian/time-of-day study; contrasts should retain ZT stratum. "
    "Consider collection ZT as a separate moderator."
)

# ------------------------------------------------------------
# GSE278920
# Zhou et al., Cell Metabolism 2025
# PMID 40858101; DOI 10.1016/j.cmet.2025.07.010
#
# Publication Methods:
# - male C57BL/6J and 129S1/SvImJ
# - DIO: 60% HFD for 12 weeks starting at 8 weeks of age
# - age-matched controls: normal chow
# ------------------------------------------------------------
gse = "GSE278920"
idx = set_pub(
    gse,
    "40858101",
    "10.1016/j.cmet.2025.07.010",
    "Genetics-nutrition interactions control diurnal enhancer-promoter dynamics and liver lipid metabolism",
)

citation = (
    "Zhou et al. Cell Metab. 2025;37:1961-1979.e7. "
    "PMID:40858101; DOI:10.1016/j.cmet.2025.07.010"
)

fill_var(
    idx, "sex", "MALE",
    "PUBLICATION_METHODS", citation,
    "Wild-type C57BL/6J and 129S1/SvImJ male mice were used.",
    "HIGH", "CONSISTENT",
)

fill_var(
    idx, "hfd_duration_weeks", "12",
    "PUBLICATION_METHODS", citation,
    "Diet-induced obesity mice were fed HFD for 12 weeks starting at 8 weeks of age.",
    "HIGH", "MANUAL_CONFIRMED",
)

fill_var(
    idx, "hfd_fat_kcal_percent", "60",
    "PUBLICATION_METHODS", citation,
    "Nutritional challenge used a 60% high-fat rodent diet.",
    "HIGH", "MANUAL_CONFIRMED",
)

fill_var(
    idx, "control_diet_type", "CHOW",
    "PUBLICATION_METHODS", citation,
    "Age-matched control mice were fed normal chow.",
    "HIGH", "CONSISTENT",
)

idx_129 = idx & x["Contrast_id"].isin(
    [f"GSE278920__C{i:02d}" for i in range(1, 7)]
)
fill_var(
    idx_129, "strain", "129S1/SvImJ",
    "PUBLICATION_METHODS", citation,
    "Selected 129 strain strata correspond to wild-type 129S1/SvImJ mice.",
    "HIGH", "MANUAL_CONFIRMED",
    final_note="Refines GEO-auto coarse normalized value '129' to exact substrain.",
)

idx_b6 = idx & x["Contrast_id"].isin(
    [f"GSE278920__C{i:02d}" for i in range(7, 13)]
)
fill_var(
    idx_b6, "strain", "C57BL/6J",
    "PUBLICATION_METHODS", citation,
    "Selected B6 strata correspond to wild-type C57BL/6J mice.",
    "HIGH", "CONSISTENT",
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "GEO design and publication explicitly compare normal chow vs HFD "
    "within C57BL/6J and 129S1/SvImJ backgrounds."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH01_ADJUDICATED"
x.loc[idx, "general_note"] = (
    "Six circadian sampling times per strain; retain strain and ZT as context. "
    "Within-GSE contrasts are dependent."
)

# ------------------------------------------------------------
# GSE236201
# Yao et al., Cell Metabolism 2025
# PMID 39814018; DOI 10.1016/j.cmet.2024.12.002
#
# GEO overall design:
# - 4-week-old male offspring
# - low-fat diet vs high-fat diet for 12 weeks
#
# Publication Methods:
# - C57BL/6J breeding colony
# - HFD D12492, 60 kcal% fat
# - LFD D12450J
# - 12 weeks
# ------------------------------------------------------------
gse = "GSE236201"
idx = set_pub(
    gse,
    "39814018",
    "10.1016/j.cmet.2024.12.002",
    "Maternal circadian rhythms during pregnancy dictate metabolic plasticity in offspring",
)

citation = (
    "Yao et al. Cell Metab. 2025;37:395-412.e6. "
    "PMID:39814018; DOI:10.1016/j.cmet.2024.12.002"
)

fill_var(
    idx, "sex", "MALE",
    "GEO_SERIES", "GSE236201 Overall design + publication Methods",
    "GEO RNA-seq overall design specifies 4-week-old male offspring; "
    "publication describes the offspring dietary challenge.",
    "HIGH", "MANUAL_CONFIRMED",
)

fill_var(
    idx, "strain", "C57BL/6J",
    "PUBLICATION_METHODS", citation,
    "Parental/experimental mice are C57BL/6J; selected offspring RNA-seq "
    "cohorts derive from this C57BL/6J design.",
    "HIGH", "MANUAL_CONFIRMED",
)

fill_var(
    idx, "hfd_duration_weeks", "12",
    "PUBLICATION_METHODS", citation,
    "Four-week-old offspring were fed HFD or LFD for 12 weeks.",
    "HIGH", "CONSISTENT",
)

fill_var(
    idx, "hfd_fat_kcal_percent", "60",
    "PUBLICATION_METHODS", citation,
    "HFD: D12492, Research Diets, 60 kcal% fat.",
    "HIGH", "MANUAL_CONFIRMED",
)

fill_var(
    idx, "control_diet_type", "LFD",
    "PUBLICATION_METHODS", citation,
    "Low-fat control diet: Research Diets D12450J.",
    "HIGH", "CONSISTENT",
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Selected GSMs explicitly identify low-fat vs high-fat diet and "
    "maternal control vs circadian-disruption strata at matched ZT."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH01_ADJUDICATED"
x.loc[idx, "general_note"] = (
    "Maternal circadian condition and ZT are preserved strata. "
    "Consider collection ZT as a separate moderator/context variable."
)

# ------------------------------------------------------------
# QA
# ------------------------------------------------------------
batch_gses = ["GSE108688", "GSE278920", "GSE236201"]
b = x[x["GSE"].isin(batch_gses)].copy()

assert b["Contrast_id"].nunique() == 30
assert b["GSE"].nunique() == 3

for var in [
    "hfd_duration_weeks",
    "sex",
    "strain",
    "hfd_fat_kcal_percent",
    "control_diet_type",
]:
    assert (b[f"{var}_final"] != "").all(), var
    assert (b[f"{var}_adjudication_status"] != "NOT_REVIEWED").all(), var

# Ensure GSE108688 parser correction took effect.
z = b[b["GSE"] == "GSE108688"]
assert set(z["hfd_fat_kcal_percent_final"]) == {"58"}
assert set(z["control_diet_type_final"]) == {"LFD"}
assert set(z["sex_final"]) == {"MALE"}
assert set(z["hfd_duration_weeks_final"]) == {"12"}

# Ensure GSE278920 strain structure.
z = b[b["GSE"] == "GSE278920"]
assert set(
    z[z["Contrast_id"].isin([f"GSE278920__C{i:02d}" for i in range(1,7)])]["strain_final"]
) == {"129S1/SvImJ"}
assert set(
    z[z["Contrast_id"].isin([f"GSE278920__C{i:02d}" for i in range(7,13)])]["strain_final"]
) == {"C57BL/6J"}

# Ensure GSE236201 final values.
z = b[b["GSE"] == "GSE236201"]
assert set(z["sex_final"]) == {"MALE"}
assert set(z["strain_final"]) == {"C57BL/6J"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"60"}
assert set(z["control_diet_type_final"]) == {"LFD"}

x.to_csv(OUT, sep="\t", index=False)

print("===== HFD MODERATOR MANUAL ADJUDICATION BATCH01 PASS =====")
print("GSE adjudicated:", b["GSE"].nunique())
print("Contrasts adjudicated:", b["Contrast_id"].nunique())

print("\n===== BATCH01 FINAL SUMMARY =====")
cols = [
    "GSE", "Contrast_id",
    "sex_final", "strain_final",
    "hfd_duration_weeks_final",
    "hfd_fat_kcal_percent_final",
    "control_diet_type_final",
]
print(b[cols].to_string(index=False))

print("\nSaved:")
print(OUT.relative_to(ROOT))
if VOCAB.exists():
    print(VOCAB_OUT.relative_to(ROOT))

print("\nNOTE:")
print("GSE108688 HFD fat 11% was an AUTO-EXTRACTION error:")
print("11% = control D12329; HFD D12331 = 58 kcal% fat.")
print("No frozen contrast eligibility was changed.")
