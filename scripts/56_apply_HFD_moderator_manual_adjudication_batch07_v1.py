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

MASTER = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch06_v1_1.tsv"
GSM_AUTO = (
    ROOT /
    "00_metadata/contrasts/moderators/"
    "HFD_primary837_GSM_moderator_extraction_GEOauto_v1_1_1.tsv"
)

OUT = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch07_v1.tsv"
ATTEMPT_OUT = OUTDIR / "HFD_batch07_source_insufficient_review_log_v1.tsv"
GSE324_AUDIT_OUT = OUTDIR / "HFD_GSE324134_selected_metadata_audit_batch07_v1.tsv"

for p in [MASTER, GSM_AUTO]:
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

def selected_df(gse, cid=None):
    z = g[g["GSE"].eq(gse)].copy()
    if cid is not None:
        z = z[z["Contrast_id"].eq(cid)].copy()
    return z

def selected_text(gse, cid=None):
    z = selected_df(gse, cid)
    vals = []
    for c in ["sample_title_raw", "sample_characteristics_raw"]:
        if c in z.columns:
            vals.extend(v for v in z[c].tolist() if v)
    return " || ".join(vals)

def infer_sex_from_text(text):
    vals = set(
        v.upper()
        for v in re.findall(
            r"(?:sex|gender)\s*:\s*(male|female)",
            text,
            flags=re.I,
        )
    )
    if vals == {"MALE"}:
        return "MALE"
    if vals == {"FEMALE"}:
        return "FEMALE"
    if vals == {"MALE", "FEMALE"}:
        return "MIXED"
    return ""

def canon_strain_from_text(text):
    # Only exact explicit tokens; never infer J/N substrain.
    if re.search(r"C57BL[/ ]?6J", text, re.I):
        return "C57BL/6J"
    if re.search(r"C57BL[/ ]?6N", text, re.I):
        return "C57BL/6N"
    if re.search(r"C57BL[/ ]?6\b", text, re.I):
        return "C57BL/6"
    return ""

attempts = []

# ============================================================
# GSE267123 — duration strict upgrade; likely COMPLETE
# Ramirez et al. JCI Insight. 2025.
# PMID 39998893
#
# GEO:
# - C57BL/6J
# - male + female
# - FFD/high-fat vs control
# - n=6 per diet per sex
# - 20 weeks
#
# Publication:
# - GAN diet D09100310 = 40% fat, 20% fructose, 2% cholesterol
# - matched GAN control D09100304 = 10% fat
#
# If all 24 Series GSM are in the one frozen contrast, sex=MIXED.
# If frozen selection is a subset and sample-specific sex is unavailable,
# sex stays unresolved rather than being inferred.
# ============================================================
gse = "GSE267123"
idx = set_pub(
    gse,
    "39998893",
    "",
    "The transcription factor ZNF469 regulates collagen production in liver fibrosis",
)
citation = (
    "Ramirez-Solano et al. JCI Insight. 2025. PMID:39998893; "
    "GEO:GSE267123"
)
z = selected_df(gse)
t = selected_text(gse)

if z["Contrast_id"].nunique() != 1:
    raise RuntimeError(
        f"GSE267123: expected one frozen contrast, observed "
        f"{z['Contrast_id'].nunique()}."
    )
if re.search(r"\bMCD\b|\bob/?ob\b|ethionine", t, re.I):
    raise RuntimeError(
        "GSE267123: non-FFD liver-fibrosis model samples detected in frozen contrast."
    )

fill_var(
    idx, "strain", "C57BL/6J",
    "GEO_SERIES+PUBLICATION", citation,
    "GEO overall design and publication identify the FFD/GAN cohort as C57BL/6J."
)
fill_var(
    idx, "hfd_duration_weeks", "20",
    "GEO_SERIES", "GSE267123 GEO overall design",
    "Male and female C57BL/6J animals were fed FFD/high-fat or control diet for 20 weeks."
)
fill_var(
    idx, "hfd_fat_kcal_percent", "40",
    "PUBLICATION_METHODS", citation,
    "The corresponding GAN diet was Research Diets D09100310, 40% fat "
    "(palm oil), 20% fructose, and 2% cholesterol."
)
fill_var(
    idx, "control_diet_type", "LFD",
    "PUBLICATION_METHODS", citation,
    "Matched GAN control was Research Diets D09100304; published diet "
    "composition identifies it as a defined 10 kcal% fat control."
)

n267 = z["GSM"].nunique()
sex267_direct = infer_sex_from_text(t)

if sex267_direct:
    sex267 = sex267_direct
    sex_evidence = (
        "Sex resolved directly from selected frozen GSM metadata."
    )
    fill_var(
        idx, "sex", sex267, "GEO_GSM", "GSE267123 selected GSM metadata",
        sex_evidence
    )
elif n267 == 24:
    # GEO Series has exactly 24 samples: n=6 per diet per sex.
    sex267 = "MIXED"
    fill_var(
        idx, "sex", "MIXED", "GEO_SERIES", "GSE267123 GEO overall design",
        "All 24 Series GSM are present in the single frozen contrast. "
        "GEO design states n=6 per diet per sex; therefore the contrast "
        "contains both male and female mice.",
        final_note=(
            "Known mixed-sex contrast; retain in overall meta-analysis but "
            "exclude from binary male-vs-female subgroup analysis."
        ),
    )
else:
    sex267 = ""
    mark_not_reported(
        idx, "sex", "GEO_SERIES", "GSE267123 GEO Series/GSM metadata",
        "GEO states both sexes were studied, but the frozen selected subset "
        "does not expose sample-specific sex sufficiently to assign this "
        "single contrast as male, female, or mixed without inference."
    )

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    f"Frozen selected GSE267123 universe contains {n267} GSM in one FFD-vs-control contrast; "
    "non-FFD fibrosis-model samples are excluded by QA."
)
x.loc[idx, "overall_adjudication_status"] = (
    "BATCH07_ADJUDICATED" if sex267 else "BATCH07_PARTIAL_REVIEWED"
)

if not sex267:
    attempts.append({
        "GSE": gse,
        "unresolved_fields": "sex",
        "reason": (
            "Series includes both sexes but selected-subset sex cannot be "
            "assigned without sample-specific evidence."
        ),
    })

# ============================================================
# GSE283060 — COMPLETE
# Corbalan et al. PLoS One. 2025.
# PMID 41091718; DOI 10.1371/journal.pone.0334213
#
# - 8-week-old male
# - C57BL/6J background
# - 16 weeks
# - chow Picochow 5053 = 13 kcal% fat
# - Western HFD RD12079B = 40 kcal% fat
#
# Frozen primary must be Group 1 WT HFD vs Group 2 WT chow.
# Humanized groups 3/4 must not be present.
# ============================================================
gse = "GSE283060"
idx = set_pub(
    gse,
    "41091718",
    "10.1371/journal.pone.0334213",
    "Humanized monoacylglycerol acyltransferase 2 mice on a high-fat diet exhibit impaired liver detoxification during metabolic dysfunction-associated steatotic liver disease",
)
citation = (
    "Corbalan et al. PLoS One. 2025. PMID:41091718; "
    "DOI:10.1371/journal.pone.0334213; GEO:GSE283060"
)
z = selected_df(gse)
t = selected_text(gse)

if z["Contrast_id"].nunique() != 1:
    raise RuntimeError("GSE283060: expected one frozen primary contrast.")
if re.search(r"Group\s*[34]\b|Humanized|HuMgat2", t, re.I):
    raise RuntimeError(
        "GSE283060: humanized Mgat2 samples detected in frozen primary contrast; "
        "expected WT Group 1 vs Group 2."
    )

fill_var(
    idx, "sex", "MALE", "PUBLICATION_METHODS", citation,
    "Metabolic feeding study used 8-week-old male mMgat2 and HuMgat2 mice."
)
fill_var(
    idx, "strain", "C57BL/6J", "PUBLICATION_METHODS", citation,
    "The mouse model is maintained/generated on a C57BL/6J background; "
    "frozen primary selected samples are WT mMgat2."
)
fill_var(
    idx, "hfd_duration_weeks", "16", "GEO_SERIES+PUBLICATION", citation,
    "WT and humanized mice were fed chow or HFD for 16 weeks; RNA-seq was performed at study end."
)
fill_var(
    idx, "hfd_fat_kcal_percent", "40", "PUBLICATION_METHODS", citation,
    "HFD was Research Diets RD12079B Western Diet, 40% calories from fat."
)
fill_var(
    idx, "control_diet_type", "CHOW", "PUBLICATION_METHODS", citation,
    "Control was Purina Picochow 5053, a standard chow diet (13% calories from fat)."
)
x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected samples must be WT Group 1 HFD vs WT Group 2 chow; "
    "humanized Groups 3/4 are rejected by QA."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH07_ADJUDICATED"
x.loc[idx, "general_note"] = "Mice were not fasted before euthanasia."

# ============================================================
# GSE274914 — COMPLETE, major age-vs-duration correction
# Morral et al. Hepatol Commun. 2024.
# PMID 39185904
#
# GEO auto duration 52 weeks is age-at-start contamination.
# Actual HFD/LFD feeding duration = 10 weeks for ALL four strata.
# - C57BL/6J
# - male/female
# - young start 7 weeks, aged start 52 weeks
# - HFD D12492 = 60 kcal% fat
# - matched LFD D12450K = 10 kcal% fat
# ============================================================
gse = "GSE274914"
idx = set_pub(
    gse,
    "39185904",
    "",
    "Sex- and age-associated factors drive the pathophysiology of MASLD",
)
citation = (
    "Morral et al. Hepatol Commun. 2024. PMID:39185904; "
    "GEO:GSE274914"
)

cids = sorted(x.loc[idx, "Contrast_id"].unique())
if len(cids) != 4:
    raise RuntimeError(
        f"GSE274914: expected four frozen contrasts, observed {len(cids)}."
    )

for cid in cids:
    cidx = x["Contrast_id"].eq(cid)
    ct = selected_text(gse, cid)
    csex = infer_sex_from_text(ct)

    if csex not in {"MALE", "FEMALE"}:
        raise RuntimeError(
            f"{cid}: could not resolve a single sex from selected GSM metadata."
        )

    fill_var(
        cidx, "sex", csex, "GEO_GSM",
        "GSE274914 selected GSM characteristics + publication",
        f"Selected frozen GSM metadata resolve this age/diet stratum as {csex.lower()}."
    )
    fill_var(
        cidx, "strain", "C57BL/6J",
        "PUBLICATION_METHODS", citation,
        "Male and female young and aged mice were C57BL/6J."
    )
    fill_var(
        cidx, "hfd_duration_weeks", "10",
        "PUBLICATION_METHODS", citation,
        "Young (7-week-old) and aged (52-week-old) mice were fed HFD or matched LFD for 10 weeks.",
        status="AUTO_EXTRACTION_ERROR_CORRECTED",
        final_note=(
            "GEO-auto 52 weeks represented age at diet initiation, not HFD exposure. "
            "All four frozen age/sex strata have 10 weeks diet exposure."
        ),
    )
    fill_var(
        cidx, "hfd_fat_kcal_percent", "60",
        "PUBLICATION_METHODS", citation,
        "HFD was Research Diets D12492, 60 kcal% fat."
    )
    fill_var(
        cidx, "control_diet_type", "LFD",
        "PUBLICATION_METHODS", citation,
        "Matched control was Research Diets D12450K, 10 kcal% fat."
    )

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Four frozen contrasts preserve age-by-sex HFD-vs-LFD strata separately; "
    "diet duration is 10 weeks in each."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH07_ADJUDICATED"
x.loc[idx, "general_note"] = (
    "Age-at-start (7 vs 52 weeks) is a context variable and must not be used as HFD duration."
)

# ============================================================
# GSE197649 — COMPLETE
# Morral et al. BBA Gene Regul Mech. 2021.
# PMID 33556624; DOI 10.1016/j.bbagrm.2021.194691
#
# - control mice from colony maintained on C57BL/6J
# - 8-week-old male
# - 10 weeks
# - HFD D12492 = 60 kcal% fat
# - control 2018SX = crude laboratory chow, 18 kcal% fat
# ============================================================
gse = "GSE197649"
idx = set_pub(
    gse,
    "33556624",
    "10.1016/j.bbagrm.2021.194691",
    "Aberrant gene expression induced by a high fat diet is linked to H3K9 acetylation in the promoter-proximal region",
)
citation = (
    "Morral et al. Biochim Biophys Acta Gene Regul Mech. 2021;1864:194691. "
    "PMID:33556624; DOI:10.1016/j.bbagrm.2021.194691"
)
z = selected_df(gse)
if z["Contrast_id"].nunique() != 1:
    raise RuntimeError("GSE197649: expected one frozen primary contrast.")

fill_var(
    idx, "sex", "MALE", "PUBLICATION_METHODS", citation,
    "Eight-week-old male control mice were used for the HFD/chow experiment."
)
fill_var(
    idx, "strain", "C57BL/6J", "PUBLICATION_METHODS", citation,
    "The control-mouse colony is maintained on a C57BL/6J background."
)
fill_var(
    idx, "hfd_duration_weeks", "10", "PUBLICATION_METHODS", citation,
    "Mice were fed control or HFD for 10 weeks."
)
fill_var(
    idx, "hfd_fat_kcal_percent", "60", "PUBLICATION_METHODS", citation,
    "HFD was Research Diets D12492, 60 kcal% fat."
)
fill_var(
    idx, "control_diet_type", "CHOW", "PUBLICATION_METHODS", citation,
    "Control was Envigo 2018SX, explicitly described as crude laboratory chow (18% kcal fat)."
)
x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen primary is the control-genotype chow-vs-HFD liver comparison."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH07_ADJUDICATED"

# ============================================================
# GSE288077 — PARTIAL, duration strict upgrade
# Han et al. J Hepatol. 2025.
# PMID 39900119
#
# GEO:
# - Control1-4 = C57BL/6, chow
# - MASLD1-4 = C57BL/6, HFD
# - HFD maintained 16 weeks
# HBV and HBV-MASLD samples must NOT enter frozen primary.
#
# Sex and kcal% fat are not explicitly reported in the accessible
# study/GEO evidence -> leave unresolved.
# ============================================================
gse = "GSE288077"
idx = set_pub(
    gse,
    "39900119",
    "",
    "Metabolic dysfunction-associated steatotic liver disease may increase intrahepatic interferon gene signatures in patients with chronic hepatitis B",
)
citation = "Han et al. J Hepatol. 2025. PMID:39900119; GEO:GSE288077"
z = selected_df(gse)
t = selected_text(gse)

if z["Contrast_id"].nunique() != 1:
    raise RuntimeError("GSE288077: expected one frozen primary contrast.")
if re.search(r"\bHBV\b|HBV[-_ ]?MASLD|rAAV", t, re.I):
    raise RuntimeError(
        "GSE288077: HBV/rAAV samples detected in frozen primary contrast."
    )
if not re.search(r"Control[1-4]", t, re.I):
    raise RuntimeError("GSE288077: Control1-4 samples not detectable.")
if not re.search(r"MASLD[1-4]", t, re.I):
    raise RuntimeError("GSE288077: MASLD1-4 samples not detectable.")

fill_var(
    idx, "strain", "C57BL/6", "GEO_GSM", "GSE288077 selected GSM characteristics",
    "Selected Control and MASLD mouse liver GSM explicitly report strain C57BL/6."
)
fill_var(
    idx, "hfd_duration_weeks", "16", "GEO_SERIES+PUBLICATION", citation,
    "GEO overall design states HFD was maintained for 16 weeks to establish the MASLD model."
)
fill_var(
    idx, "control_diet_type", "CHOW", "GEO_GSM", "GSE288077 selected Control GSM",
    "Selected Control1-4 mouse liver samples are explicitly labeled chow diet."
)
mark_not_reported(
    idx, "sex", "GEO_SERIES+PUBLICATION", citation,
    "Sex of the selected Control/MASLD RNA-seq cohort was not explicitly reported in the accessible GEO/publication evidence."
)
mark_not_reported(
    idx, "hfd_fat_kcal_percent", "GEO_SERIES+PUBLICATION", citation,
    "The accessible GEO/publication evidence identifies HFD but does not explicitly report HFD kcal% fat for the selected Control/MASLD cohort."
)

x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen primary is Control1-4 vs MASLD1-4; HBV and HBV-MASLD cohorts are rejected by QA."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH07_PARTIAL_REVIEWED"

attempts.append({
    "GSE": gse,
    "unresolved_fields": "sex,hfd_fat_kcal_percent",
    "reason": (
        "Accessible GEO/publication evidence does not explicitly report these "
        "fields for the selected Control/MASLD mouse RNA-seq cohort."
    ),
})

# ============================================================
# GSE324134 — CONDITIONAL DIRECT-GEO adjudication only
#
# Public indexing is incomplete. Do NOT assume that the related HQHF
# publications describe this exact deposited RNA-seq cohort.
#
# We inspect the frozen selected GSM metadata:
# - if duration=12 weeks is directly present, upgrade duration to HIGH;
# - sex/strain/control are filled only if directly explicit;
# - ingredient wt% is NOT converted to kcal% fat;
# - if direct evidence is insufficient, leave this GSE untouched so the
#   existing GEO_AUTO_LOW record remains auditable and re-reviewable.
# ============================================================
gse = "GSE324134"
z324 = selected_df(gse)
t324 = selected_text(gse)

audit_cols = [
    c for c in [
        "GSE", "Contrast_id", "GSM",
        "sample_title_raw", "sample_characteristics_raw"
    ] if c in z324.columns
]
z324[audit_cols].to_csv(GSE324_AUDIT_OUT, sep="\t", index=False)

print("\n===== GSE324134 SELECTED METADATA AUDIT =====")
print("selected_GSM:", z324["GSM"].nunique())
print("selected_contrasts:", z324["Contrast_id"].nunique())
print("selected_text_preview:", t324[:2500])

duration324_direct = bool(
    re.search(r"\b12\s*(?:weeks?|wks?|wk)\b", t324, re.I)
)
sex324 = infer_sex_from_text(t324)
strain324 = canon_strain_from_text(t324)
has_control324 = bool(
    re.search(r"regular\s*(?:feed|diet)|normal\s*(?:feed|diet)|chow", t324, re.I)
)
has_hfd324 = bool(re.search(r"high[- ]?fat|\bHFD\b", t324, re.I))
has_hqhf324 = bool(re.search(r"\bHQHF\b|Huatan", t324, re.I))

gse324_adjudicated = False

if has_hqhf324:
    raise RuntimeError(
        "GSE324134: HQHF intervention samples detected in frozen primary contrast."
    )

if duration324_direct and has_hfd324:
    idx = x["GSE"].eq(gse)

    fill_var(
        idx, "hfd_duration_weeks", "12",
        "GEO_GSM", "GSE324134 selected frozen GSM metadata",
        "Selected frozen GSM metadata directly report 12 weeks of HFD exposure."
    )

    if sex324:
        fill_var(
            idx, "sex", sex324,
            "GEO_GSM", "GSE324134 selected frozen GSM metadata",
            "Sex resolved directly from selected frozen GSM metadata."
        )
    else:
        mark_not_reported(
            idx, "sex", "GEO_GSM", "GSE324134 selected frozen GSM metadata",
            "Selected frozen GSM metadata do not explicitly resolve sex."
        )

    if strain324:
        fill_var(
            idx, "strain", strain324,
            "GEO_GSM", "GSE324134 selected frozen GSM metadata",
            "Strain resolved directly from selected frozen GSM metadata."
        )
    else:
        mark_not_reported(
            idx, "strain", "GEO_GSM", "GSE324134 selected frozen GSM metadata",
            "Selected frozen GSM metadata do not explicitly resolve strain."
        )

    if has_control324:
        fill_var(
            idx, "control_diet_type", "CHOW",
            "GEO_GSM", "GSE324134 selected frozen GSM metadata",
            "Selected control samples explicitly identify regular/normal feed or chow."
        )
    else:
        fill_var(
            idx, "control_diet_type", "CONTROL_DIET_UNSPECIFIED",
            "GEO_GSM", "GSE324134 selected frozen GSM metadata",
            "A clean control-vs-HFD contrast is frozen, but control diet composition/type "
            "is not explicit enough to classify as chow versus purified LFD."
        )

    # Never calculate kcal% fat from ingredient weight percentages.
    mark_not_reported(
        idx, "hfd_fat_kcal_percent",
        "GEO_GSM", "GSE324134 selected frozen GSM metadata",
        "HFD may be reported by ingredient weight composition; no explicit kcal% fat "
        "was available and no conversion from ingredient wt% to kcal% was attempted."
    )

    x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
    x.loc[idx, "sample_role_identity_evidence"] = (
        "Frozen selected metadata directly support HFD duration; HQHF intervention "
        "samples are absent by QA."
    )
    x.loc[idx, "overall_adjudication_status"] = "BATCH07_PARTIAL_REVIEWED"
    gse324_adjudicated = True

    unresolved324 = []
    if not sex324:
        unresolved324.append("sex")
    if not strain324:
        unresolved324.append("strain")
    unresolved324.append("hfd_fat_kcal_percent")

    attempts.append({
        "GSE": gse,
        "unresolved_fields": ",".join(unresolved324),
        "reason": (
            "Only directly supported frozen-GSM metadata were adjudicated; "
            "ingredient weight percentages were not converted to kcal% fat."
        ),
    })
else:
    # Intentionally DO NOT modify master fields/status. Existing LOW auto
    # evidence remains preserved for later source review.
    attempts.append({
        "GSE": gse,
        "unresolved_fields": "BATCH07_DEFERRED",
        "reason": (
            "Selected frozen GSM metadata did not directly contain sufficient "
            "12-week HFD evidence to supersede GEO_AUTO_LOW. No related HQHF "
            "publication was assumed to describe this exact deposited cohort."
        ),
    })

# ------------------------------------------------------------
# QA
# ------------------------------------------------------------
assert x["GSE"].nunique() == 64
assert x["Contrast_id"].nunique() == 105

# Definite strict-duration upgrades
for gse, dur in {
    "GSE267123": "20",
    "GSE283060": "16",
    "GSE274914": "10",
    "GSE197649": "10",
    "GSE288077": "16",
}.items():
    z = x[x["GSE"].eq(gse)]
    assert set(z["hfd_duration_weeks_final"]) == {dur}, (gse, dur)
    assert set(z["hfd_duration_weeks_final_confidence"]) == {"HIGH"}, gse

# GSE274914 all 4 are 10-week strata.
z = x[x["GSE"].eq("GSE274914")]
assert z["Contrast_id"].nunique() == 4
assert set(z["hfd_duration_weeks_final"]) == {"10"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"60"}
assert set(z["control_diet_type_final"]) == {"LFD"}
assert set(z["strain_final"]) == {"C57BL/6J"}
assert set(z["sex_final"]).issubset({"MALE", "FEMALE"})
assert len(set(z["sex_final"])) == 2

# GSE197649 exact
z = x[x["GSE"].eq("GSE197649")]
assert set(z["sex_final"]) == {"MALE"}
assert set(z["strain_final"]) == {"C57BL/6J"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"60"}
assert set(z["control_diet_type_final"]) == {"CHOW"}

# GSE283060 exact
z = x[x["GSE"].eq("GSE283060")]
assert set(z["sex_final"]) == {"MALE"}
assert set(z["strain_final"]) == {"C57BL/6J"}
assert set(z["hfd_fat_kcal_percent_final"]) == {"40"}
assert set(z["control_diet_type_final"]) == {"CHOW"}

# GSE288077 source-limited
z = x[x["GSE"].eq("GSE288077")]
assert set(z["strain_final"]) == {"C57BL/6"}
assert set(z["sex_final"]) == {""}
assert set(z["hfd_fat_kcal_percent_final"]) == {""}
assert set(z["control_diet_type_final"]) == {"CHOW"}

# Preserve all previously reviewed statuses exactly.
prev = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
old_reviewed = prev[
    prev["overall_adjudication_status"].str.contains(
        "ADJUDICATED|PARTIAL_REVIEWED", regex=True, na=False
    )
]["GSE"].unique()

for old_gse in old_reviewed:
    a = prev[prev["GSE"].eq(old_gse)].sort_values("Contrast_id")
    bb = x[x["GSE"].eq(old_gse)].sort_values("Contrast_id")
    assert (
        a["overall_adjudication_status"].tolist()
        == bb["overall_adjudication_status"].tolist()
    ), old_gse

x.to_csv(OUT, sep="\t", index=False)
pd.DataFrame(attempts).to_csv(ATTEMPT_OUT, sep="\t", index=False)

# Status summary.
status_rows = []
for gg, zz in x.groupby("GSE"):
    sts = set(s for s in zz["overall_adjudication_status"] if s)
    if any("PARTIAL_REVIEWED" in s for s in sts):
        st = "PARTIAL_REVIEWED"
    elif sts and all("ADJUDICATED" in s for s in sts):
        st = "FULLY_ADJUDICATED"
    else:
        st = "NOT_YET_REVIEWED"
    status_rows.append((gg, st, zz["Contrast_id"].nunique()))

sdf = pd.DataFrame(status_rows, columns=["GSE", "status", "n_contrasts"])

print("\n===== HFD MODERATOR MANUAL ADJUDICATION BATCH07 PASS =====")
print("Definite duration strict-upgrade GSE:", 5)
print("GSE324134 direct-GEO adjudicated:", int(gse324_adjudicated))
print("Cumulative fully adjudicated GSE:",
      int((sdf["status"] == "FULLY_ADJUDICATED").sum()))
print("Cumulative partial-reviewed GSE:",
      int((sdf["status"] == "PARTIAL_REVIEWED").sum()))
print("Cumulative not-yet-reviewed GSE:",
      int((sdf["status"] == "NOT_YET_REVIEWED").sum()))

print("\n===== BATCH07 FINAL STATUS =====")
batch_gses = [
    "GSE267123", "GSE283060", "GSE274914",
    "GSE197649", "GSE288077", "GSE324134",
]
cols = [
    "GSE", "Contrast_id", "overall_adjudication_status",
    "sex_final", "strain_final", "hfd_duration_weeks_final",
    "hfd_duration_weeks_final_confidence",
    "hfd_fat_kcal_percent_final", "control_diet_type_final",
]
print(
    x[x["GSE"].isin(batch_gses)][cols]
    .sort_values(["GSE", "Contrast_id"])
    .to_string(index=False)
)

print("\n===== SOURCE-INSUFFICIENT / DEFERRED =====")
print(pd.DataFrame(attempts).to_string(index=False))

print("\nKEY CORRECTIONS:")
print("- GSE274914: auto duration=52 was age at diet initiation;")
print("  all 4 frozen sex/age strata have 10 weeks of HFD/LFD feeding.")
print("- GSE267123: FFD/GAN diet = 40% fat; matched D09100304 control = LFD.")
print("- GSE288077: duration=16 and C57BL/6/chow are confirmed; sex/fat remain unknown.")
print("- GSE324134: no related publication is assumed to describe the deposited cohort;")
print("  it is upgraded only if selected frozen GSM metadata directly support 12 weeks.")

print("\nEXPECTED DURATION EFFECT:")
print("- Definite strict duration independent GSE: 28 -> 33.")
if gse324_adjudicated:
    print("- GSE324134 also upgraded directly from frozen GSM evidence: expected 34.")
else:
    print("- GSE324134 deferred: expected strict duration remains 33.")

print("\nSaved:")
print(OUT.relative_to(ROOT))
print(ATTEMPT_OUT.relative_to(ROOT))
print(GSE324_AUDIT_OUT.relative_to(ROOT))

print("\nIMPORTANT:")
print("- Ingredient weight percentages are never converted to kcal% fat.")
print("- Age-at-start is not duration.")
print("- Mixed sex is retained as MIXED, never majority-labeled.")
print("- No frozen eligibility or contrast membership changed.")
