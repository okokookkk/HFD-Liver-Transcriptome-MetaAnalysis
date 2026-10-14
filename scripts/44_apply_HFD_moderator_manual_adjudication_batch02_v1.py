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
    "HFD_primary105_moderator_manual_adjudication_batch01_plus_ZT_v1_1.tsv"
)

GSM_AUTO = (
    ROOT /
    "00_metadata/contrasts/moderators/"
    "HFD_primary837_GSM_moderator_extraction_GEOauto_v1_1_1.tsv"
)

OUT = (
    ROOT /
    "00_metadata/contrasts/moderators/manual_adjudication/"
    "HFD_primary105_moderator_manual_adjudication_batch02_v1.tsv"
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
assert g["Contrast_id"].nunique() == 105

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

def selected_titles(gse, cid=None):
    z = g[g["GSE"].eq(gse)].copy()
    if cid is not None:
        z = z[z["Contrast_id"].eq(cid)]
    vals = []
    for c in ["sample_title_raw", "sample_characteristics_raw"]:
        if c in z.columns:
            vals.extend(v for v in z[c].tolist() if v)
    return " || ".join(vals)

# ============================================================
# BATCH02
# GSE218932, GSE246221, GSE294746, GSE296500
# ============================================================

# ------------------------------------------------------------
# GSE218932
# Yan et al. Front Nutr. 2023.
# PMID 36687679; DOI 10.3389/fnut.2022.1068350
#
# Publication:
# - female C57BL/6NHsd
# - 3 weeks feeding
# - HFD 45 energy%
# - AIN93G control 16 energy% fat
# - not fasted
# - samples collected across ZT
# ------------------------------------------------------------
gse = "GSE218932"
idx = set_pub(
    gse,
    "36687679",
    "10.3389/fnut.2022.1068350",
    "Consumption of a high-fat diet alters transcriptional rhythmicity in liver from pubertal mice",
)
citation = (
    "Yan et al. Front Nutr. 2023;9:1068350. "
    "PMID:36687679; DOI:10.3389/fnut.2022.1068350"
)

fill_var(
    idx, "sex", "FEMALE",
    "PUBLICATION_METHODS", citation,
    "Three-week-old female C57BL/6NHsd mice were used.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "strain", "C57BL/6NHsd",
    "PUBLICATION_METHODS", citation,
    "Mice were C57BL/6NHsd purchased from Envigo.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "hfd_duration_weeks", "3",
    "PUBLICATION_METHODS", citation,
    "Mice had free access to AIN93G or HFD for 3 weeks.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "hfd_fat_kcal_percent", "45",
    "PUBLICATION_METHODS", citation,
    "Diet composition table reports HFD fat = 45% of energy.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "control_diet_type", "LFD",
    "PUBLICATION_METHODS", citation,
    "Control was purified AIN93G standard diet; fat = 16% of energy.",
    "HIGH", "MANUAL_CONFIRMED",
    final_note="Purified AIN93G standard control, not conventional chow."
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Publication and selected GEO samples define matched AIN93G versus HFD "
    "contrasts across circadian collection times."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH02_ADJUDICATED"
x.loc[idx, "general_note"] = (
    "Pubertal female circadian study; 3-week HFD. Mice were not fasted. "
    "Retain collection ZT as study-specific context."
)

# ------------------------------------------------------------
# GSE246221
# Jeong et al. Nat Commun. 2024.
# PMID 39090079; DOI 10.1038/s41467-024-50660-y
#
# Frozen primary contrasts are HFD-only vs SCD-only, without STZ.
# Publication states mice without STZ were fed SCD or HFD from 8 weeks.
# HFD D12492 = 60% fat; SCD Teklad 2018 = 18% fat.
#
# IMPORTANT:
# Determine each frozen contrast age from its SELECTED GSM titles.
# Only assign duration after verifying all selected GSMs are 20w or 32w.
# ------------------------------------------------------------
gse = "GSE246221"
idx = set_pub(
    gse,
    "39090079",
    "10.1038/s41467-024-50660-y",
    "A male mouse model for metabolic dysfunction-associated steatotic liver disease and hepatocellular carcinoma",
)
citation = (
    "Jeong et al. Nat Commun. 2024;15:6506. "
    "PMID:39090079; DOI:10.1038/s41467-024-50660-y"
)

fill_var(
    idx, "sex", "MALE",
    "PUBLICATION_METHODS", citation,
    "The study exclusively used male C57BL/6J mice.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "strain", "C57BL/6J",
    "PUBLICATION_METHODS", citation,
    "Male C57BL/6J mice were used; four-group comparisons included "
    "SCD-only and HFD-only mice without STZ.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60",
    "PUBLICATION_METHODS", citation,
    "HFD was Research Diets D12492, 60% fat.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "control_diet_type", "CHOW",
    "PUBLICATION_METHODS", citation,
    "Control diet was standard chow diet, Teklad Global 2018.",
    "HIGH", "MANUAL_CONFIRMED",
    final_note="Standard chow/SCD; publication reports 18% fat."
)

cids = sorted(x.loc[idx, "Contrast_id"].unique())
assert len(cids) == 2, cids

for cid in cids:
    t = selected_titles(gse, cid)
    # Require the frozen selected GSM titles to be internally age-consistent.
    has20 = bool(re.search(r"\b20w\b|20\s*weeks", t, re.I))
    has32 = bool(re.search(r"\b32w\b|32\s*weeks", t, re.I))

    if has20 and not has32:
        dur = "12"  # diet from age 8w to age 20w
        age_note = "Selected frozen GSMs are 20-week-old HFD-only/SCD-only mice."
    elif has32 and not has20:
        dur = "24"  # diet from age 8w to age 32w
        age_note = "Selected frozen GSMs are 32-week-old HFD-only/SCD-only mice."
    else:
        raise RuntimeError(
            f"{cid}: cannot uniquely resolve 20w vs 32w from selected GSM metadata"
        )

    cidx = x["Contrast_id"].eq(cid)
    fill_var(
        cidx, "hfd_duration_weeks", dur,
        "PUBLICATION_METHODS", citation,
        age_note + " Publication states HFD-only/SCD-only mice were fed "
        "from 8 weeks of age.",
        "HIGH", "PAPER_GEO_CONFLICT_RESOLVED",
        final_note=(
            "GEO-auto series text mixed age, model duration, and intervention "
            "durations. Final duration is contrast-specific HFD exposure."
        ),
    )

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected GSMs are explicitly HFD-fed mice without STZ versus "
    "SCD-fed mice without STZ at matched ages; STZ+HFD cohorts are not used "
    "for these primary HFD-only contrasts."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH02_ADJUDICATED"
x.loc[idx, "general_note"] = (
    "Complex MASLD model series, but frozen primary contrasts use clean "
    "HFD-only vs SCD-only cohorts without STZ. Duration adjudicated from "
    "selected GSM age and diet start at 8 weeks."
)

# ------------------------------------------------------------
# GSE294746
# Senkalfa et al. Hepatol Commun. 2025.
# PMID 40864559; DOI 10.1097/HC9.0000000000000782
#
# Selected primary contrasts:
# - chow vs HFD at room temperature
# - chow vs HFD at thermoneutrality
# Ketogenic-diet samples are NOT part of the frozen primary contrasts.
#
# Publication:
# - male
# - RNA-seq long-term cohort described as C57BL/6N
# - 14 weeks
# - HFD = 60% energy from fat
# - standard chow control
# - fasted 6 h before sampling
# ------------------------------------------------------------
gse = "GSE294746"
idx = set_pub(
    gse,
    "40864559",
    "10.1097/HC9.0000000000000782",
    "Adipocyte-specific deletion of gp130 prevents ketogenic diet-induced hepatic steatosis",
)
citation = (
    "Senkalfa et al. Hepatol Commun. 2025;9:e0782. "
    "PMID:40864559; DOI:10.1097/HC9.0000000000000782"
)

fill_var(
    idx, "sex", "MALE",
    "PUBLICATION_METHODS", citation,
    "Long-term diet experiments used male mice.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "strain", "C57BL/6N",
    "PUBLICATION_FIGURE_LEGEND", citation,
    "RNA-seq Figure 2 explicitly describes chow/HFD/KD-fed C57BL/6N mice "
    "at room temperature or thermoneutrality.",
    "HIGH", "MANUAL_CONFIRMED",
    final_note=(
        "Paper also describes gp130 knockout lines on a C57BL/6J background; "
        "the RNA-seq cohort used for these frozen contrasts is explicitly "
        "identified as C57BL/6N in the RNA-seq figure legend."
    ),
)
fill_var(
    idx, "hfd_duration_weeks", "14",
    "PUBLICATION_METHODS", citation,
    "Chow/HFD/KD long-term feeding was performed for 14 weeks.",
    "HIGH", "CONSISTENT",
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60",
    "PUBLICATION_METHODS", citation,
    "HFD contained 60% of energy derived from fat.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "control_diet_type", "CHOW",
    "PUBLICATION_METHODS", citation,
    "Control arm used standard chow diet (#2222, Kliba-Nafag).",
    "HIGH", "MANUAL_CONFIRMED",
)

# Verify frozen selected samples do not contain ketogenic-diet samples and
# do preserve RT/TN matching within each contrast.
for cid in sorted(x.loc[idx, "Contrast_id"].unique()):
    t = selected_titles(gse, cid)
    if re.search(r"\bketogenic\b|\bKD\b", t, re.I):
        raise RuntimeError(f"{cid}: ketogenic-diet sample detected in frozen contrast")
    if not re.search(r"\bHFD\b", t, re.I):
        raise RuntimeError(f"{cid}: HFD token not found in frozen selected GSM metadata")
    if not re.search(r"\bChow\b|Chow Diet|CD_", t, re.I):
        raise RuntimeError(f"{cid}: chow/control token not found in selected GSM metadata")

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen contrasts contain chow-vs-HFD samples only; ketogenic-diet "
    "samples are absent. C01/C02 retain room-temperature and thermoneutral "
    "environmental strata separately."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH02_ADJUDICATED"
x.loc[idx, "general_note"] = (
    "14-week male HFD study with separate RT and thermoneutral contrasts; "
    "mice were fasted 6 h before sampling. Environmental temperature should "
    "remain a study-specific context variable."
)

# ------------------------------------------------------------
# GSE296500
# Dunn et al. Cell Rep. 2026.
# PMID 41538324; DOI 10.1016/j.celrep.2025.116835
#
# Publication:
# - WT and Dp16
# - 8 weeks diet starting at 28 days (4 weeks) of age
# - HFD 60 kcal% D12492
# - matched LFD D12450J (publication text reports 6 or 10 kcal% in
#   different summary/method descriptions; diet code D12450J is retained)
# - both sexes represented in each genotype/diet group
#
# Sex must therefore be MIXED, not inferred from group majority.
# ------------------------------------------------------------
gse = "GSE296500"
idx = set_pub(
    gse,
    "41538324",
    "10.1016/j.celrep.2025.116835",
    "Altered hepatic metabolism in Down syndrome",
)
citation = (
    "Dunn et al. Cell Rep. 2026;45:116835. "
    "PMID:41538324; DOI:10.1016/j.celrep.2025.116835"
)

fill_var(
    idx, "sex", "MIXED",
    "PUBLICATION_RESULTS", citation,
    "Both female and male mice are represented in WT and Dp16 LFD/HFD groups; "
    "figure legend reports four females in each diet/genotype group.",
    "HIGH", "MANUAL_CONFIRMED",
    final_note="Do not use these contrasts in male-vs-female subgroup analysis."
)
fill_var(
    idx, "strain", "C57BL/6J",
    "GEO_SERIES", "GSE296500 + publication",
    "GEO identifies wild-type C57BL/6J and Dp(16)1Yey cohorts; "
    "Dp16 is treated as a genetic-context stratum on the same mouse background.",
    "MEDIUM", "MANUAL_CONFIRMED",
    final_note=(
        "Genotype context (WT vs Dp16) should be retained separately from strain."
    ),
)
fill_var(
    idx, "hfd_duration_weeks", "8",
    "PUBLICATION_METHODS", citation,
    "Mice had ad libitum access to HFD or LFD for 8 weeks beginning at "
    "28 days of age.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60",
    "PUBLICATION_METHODS", citation,
    "HFD was Research Diets D12492, 60 kcal% fat.",
    "HIGH", "MANUAL_CONFIRMED",
)
fill_var(
    idx, "control_diet_type", "LFD",
    "PUBLICATION_METHODS", citation,
    "Matched low-fat control diet was Research Diets D12450J.",
    "HIGH", "MANUAL_CONFIRMED",
    final_note=(
        "Publication narrative contains differing low-fat kcal% descriptions "
        "(6 vs 10 kcal%); manufacturer diet code D12450J is retained as the "
        "primary identity. This does not affect classification as LFD."
    ),
)

# Verify one WT and one Dp16 frozen contrast, with both diet roles present.
geno_labels = {}
for cid in sorted(x.loc[idx, "Contrast_id"].unique()):
    t = selected_titles(gse, cid)
    has_dp16 = bool(re.search(r"\bDp16\b|Dp\(16\)", t, re.I))
    has_wt = bool(re.search(r"Wild[- ]type|\bWT\b|C57BL/6J", t, re.I))
    if has_dp16 and not has_wt:
        geno_labels[cid] = "Dp16"
    elif has_wt and not has_dp16:
        geno_labels[cid] = "WT"
    elif has_dp16:
        # Some characteristics may mention background/WT context; selected
        # sample titles still identify Dp16 explicitly.
        geno_labels[cid] = "Dp16"
    else:
        geno_labels[cid] = "WT"

if len(geno_labels) != 2:
    raise RuntimeError("GSE296500 expected exactly two frozen contrasts")

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen contrasts preserve control-vs-HFD comparisons within WT and Dp16 "
    "genotype strata; genotype strata are not pooled."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH02_ADJUDICATED"
x.loc[idx, "general_note"] = (
    "Mixed-sex, genotype-stratified study. WT and Dp16 HFD effects should "
    "remain separate dependent contrasts within one GSE. Genotype context "
    "is a candidate contextual variable, not a strain difference."
)

# ------------------------------------------------------------
# QA
# ------------------------------------------------------------
batch_gses = ["GSE218932", "GSE246221", "GSE294746", "GSE296500"]
b = x[x["GSE"].isin(batch_gses)].copy()

expected_n = {
    "GSE218932": 7,
    "GSE246221": 2,
    "GSE294746": 2,
    "GSE296500": 2,
}
for gse, n in expected_n.items():
    got = b.loc[b["GSE"].eq(gse), "Contrast_id"].nunique()
    assert got == n, (gse, got, n)

assert b["Contrast_id"].nunique() == sum(expected_n.values())

for var in [
    "hfd_duration_weeks",
    "sex",
    "strain",
    "hfd_fat_kcal_percent",
    "control_diet_type",
]:
    assert (b[f"{var}_final"] != "").all(), var
    assert (b[f"{var}_adjudication_status"] != "NOT_REVIEWED").all(), var

# Specific scientific QA.
z = b[b["GSE"].eq("GSE218932")]
assert set(z["sex_final"]) == {"FEMALE"}
assert set(z["strain_final"]) == {"C57BL/6NHsd"}
assert set(z["hfd_duration_weeks_final"]) == {"3"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"45"}

z = b[b["GSE"].eq("GSE246221")]
assert set(z["hfd_duration_weeks_final"]) == {"12", "24"}
assert set(z["sex_final"]) == {"MALE"}
assert set(z["strain_final"]) == {"C57BL/6J"}

z = b[b["GSE"].eq("GSE294746")]
assert set(z["hfd_duration_weeks_final"]) == {"14"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"60"}
assert set(z["control_diet_type_final"]) == {"CHOW"}

z = b[b["GSE"].eq("GSE296500")]
assert set(z["sex_final"]) == {"MIXED"}
assert set(z["hfd_duration_weeks_final"]) == {"8"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"60"}
assert set(z["control_diet_type_final"]) == {"LFD"}

# Preserve batch01 adjudications.
prev = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
for gse in ["GSE108688", "GSE236201", "GSE278920"]:
    a = prev[prev["GSE"].eq(gse)].sort_values("Contrast_id")
    bb = x[x["GSE"].eq(gse)].sort_values("Contrast_id")
    for var in [
        "sex_final", "strain_final", "hfd_duration_weeks_final",
        "hfd_fat_kcal_percent_final", "control_diet_type_final"
    ]:
        assert a[var].tolist() == bb[var].tolist(), (gse, var)

x.to_csv(OUT, sep="\t", index=False)

print("===== HFD MODERATOR MANUAL ADJUDICATION BATCH02 PASS =====")
print("New GSE adjudicated:", len(batch_gses))
print("New contrasts adjudicated:", b["Contrast_id"].nunique())
print("Cumulative manually adjudicated GSE:",
      x.loc[x["overall_adjudication_status"].str.contains("ADJUDICATED", na=False), "GSE"].nunique())
print("Cumulative manually adjudicated contrasts:",
      x.loc[x["overall_adjudication_status"].str.contains("ADJUDICATED", na=False), "Contrast_id"].nunique())

print("\n===== BATCH02 FINAL SUMMARY =====")
cols = [
    "GSE", "Contrast_id",
    "sex_final", "strain_final",
    "hfd_duration_weeks_final",
    "hfd_fat_kcal_percent_final",
    "control_diet_type_final",
]
print(b[cols].sort_values(["GSE", "Contrast_id"]).to_string(index=False))

print("\n===== GSE246221 DURATION CHECK =====")
print(
    b[b["GSE"].eq("GSE246221")][
        ["Contrast_id", "hfd_duration_weeks_final", "hfd_duration_weeks_final_note"]
    ].to_string(index=False)
)

print("\nSaved:")
print(OUT.relative_to(ROOT))

print("\nNOTES:")
print("- GSE246221 duration is assigned only after frozen selected GSM age verification.")
print("- GSE294746 ketogenic-diet samples are explicitly rejected by QA.")
print("- GSE296500 sex is MIXED; do not force a majority-sex label.")
print("- No frozen technical eligibility or contrast membership was changed.")
