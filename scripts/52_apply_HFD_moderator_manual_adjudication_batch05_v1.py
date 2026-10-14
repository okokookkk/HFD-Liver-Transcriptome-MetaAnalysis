#!/usr/bin/env python3
from pathlib import Path
import argparse, re
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=__import__("os").environ.get("HFD_PROJECT_ROOT", "."))
args = ap.parse_args()

ROOT = Path(args.root)
MASTER = ROOT / "00_metadata/contrasts/moderators/manual_adjudication/HFD_primary105_moderator_manual_adjudication_batch04_v1.tsv"
GSM_AUTO = ROOT / "00_metadata/contrasts/moderators/HFD_primary837_GSM_moderator_extraction_GEOauto_v1_1_1.tsv"
OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"
OUT = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch05_v1.tsv"
ATTEMPT_OUT = OUTDIR / "HFD_batch05_source_insufficient_review_log_v1.tsv"

if not MASTER.exists():
    raise FileNotFoundError(MASTER)
if not GSM_AUTO.exists():
    raise FileNotFoundError(GSM_AUTO)

x = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
g = pd.read_csv(GSM_AUTO, sep="\t", dtype=str).fillna("")

assert x["GSE"].nunique() == 64
assert x["Contrast_id"].nunique() == 105
assert g["GSM"].nunique() == 837

def fill_var(idx, var, manual_value, source_type, citation, evidence,
             confidence="HIGH", status="MANUAL_CONFIRMED", final_note=""):
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

# GSE247985 — COMPLETE
gse = "GSE247985"
idx = set_pub(
    gse, "38720938", "10.3389/fendo.2023.1332959",
    "Maternal undernutrition results in transcript changes in male offspring that may promote resistance to high fat diet induced weight gain",
)
citation = (
    "Miles et al. Front Endocrinol (Lausanne). 2024;14:1332959. "
    "PMID:38720938; DOI:10.3389/fendo.2023.1332959"
)
t = selected_text(gse)
if not re.search(r"_LIV\b|liver", t, re.I):
    raise RuntimeError("GSE247985: selected frozen samples are not clearly liver.")
if not re.search(r"HFD", t, re.I):
    raise RuntimeError("GSE247985: HFD samples not detected.")
if not re.search(r"CTRL", t, re.I):
    raise RuntimeError("GSE247985: CTRL samples not detected.")
has_fed = bool(re.search(r"FED", t, re.I))
has_cr20 = bool(re.search(r"CR20", t, re.I))
if has_fed and has_cr20:
    raise RuntimeError(
        "GSE247985: frozen contrast appears to mix FED and CR20 maternal-history strata."
    )
if not (has_fed or has_cr20):
    raise RuntimeError("GSE247985: maternal-history stratum (FED or CR20) not detectable.")

fill_var(idx, "sex", "MALE", "PUBLICATION_RESULTS", citation,
         "Liver RNA-seq was performed on adult male offspring treated for 10 weeks.")
fill_var(idx, "strain", "FVB.129P2", "PUBLICATION_METHODS", citation,
         "The colony used sighted FVB.129P2-Pde6b+ Tyrc-ch/AntJ mice; adult offspring derive from this breeding design.")
fill_var(idx, "hfd_duration_weeks", "10", "PUBLICATION_RESULTS", citation,
         "The liver RNA-seq dataset is explicitly described as 10-week-treated male offspring.",
         final_note="A separate 16-week extension was phenotyping, not the liver RNA-seq cohort.")
fill_var(idx, "hfd_fat_kcal_percent", "45", "PUBLICATION_METHODS", citation,
         "Adult offspring received Envigo TD.06415, a 45% high-fat diet.")
fill_var(idx, "control_diet_type", "LFD", "PUBLICATION_METHODS", citation,
         "Matched control was Envigo TD.110675, 10% fat.",
         final_note="Purified 10% fat control diet; classified as LFD.")
x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected liver samples contain matched HFD and CTRL within a single maternal-history stratum."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH05_ADJUDICATED"

# GSE137365 — COMPLETE
gse = "GSE137365"
idx = set_pub(
    gse, "31845524", "10.1002/mnfr.201900982",
    "Tracing the Anti-Inflammatory Mechanism/Triggers of d-Allulose: A Profile Study of Microbiome Composition and mRNA Expression in Diet-Induced Obese Mice",
)
citation = (
    "Han et al. Mol Nutr Food Res. 2020;64:e1900982. "
    "PMID:31845524; DOI:10.1002/mnfr.201900982"
)
t = selected_text(gse)
if not re.search(r"Liver", t, re.I):
    raise RuntimeError("GSE137365: liver samples not detected.")
if not re.search(r"High[- ]fat diet", t, re.I):
    raise RuntimeError("GSE137365: HFD samples not detected.")
if not re.search(r"Normal diet", t, re.I):
    raise RuntimeError("GSE137365: normal-diet control not detected.")
if re.search(r"D[- ]?allulose|ALLU", t, re.I):
    raise RuntimeError("GSE137365: D-allulose-treated HFD samples detected in frozen primary contrast.")

fill_var(idx, "sex", "MALE", "PUBLICATION_METHODS", citation,
         "The GSE137365 mouse cohort comprised male C57BL/6J mice.")
fill_var(idx, "strain", "C57BL/6J", "PUBLICATION_METHODS", citation,
         "The deposited GSE137365 cohort is explicitly C57BL/6J.",
         status="GEO_TYPO_CORRECTED",
         final_note="Corrects typo-like GEO-auto C56BL/6J using publication-supported C57BL/6J.")
fill_var(idx, "hfd_duration_weeks", "16", "PUBLICATION_METHODS", citation,
         "Mice were fed ND, HFD, or HFD+D-allulose for 16 weeks.")
fill_var(idx, "hfd_fat_kcal_percent", "39.5", "PUBLICATION_METHODS", citation,
         "The HFD provided 39.5% of total energy from fat.",
         status="AUTO_EXTRACTION_ERROR_CORRECTED",
         final_note="GEO-auto 5% was D-allulose supplementation (w/w), not HFD fat.")
fill_var(idx, "control_diet_type", "NORMAL_DIET", "PUBLICATION_METHODS", citation,
         "Control was an ND based on the AIN-76 semisynthetic diet.",
         final_note="Retained as NORMAL_DIET rather than forcing CHOW or LFD.")
x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected samples are Liver Normal diet versus Liver High-fat diet; D-allulose arm excluded."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH05_ADJUDICATED"

# GSE207281 — PARTIAL
gse = "GSE207281"
idx = set_pub(
    gse, "37138676", "10.1016/j.jhepr.2023.100726",
    "Inhibition of phospholipase D1 ameliorates hepatocyte steatosis and non-alcoholic fatty liver disease",
)
citation = (
    "Wang et al. JHEP Rep. 2023;5:100726. "
    "PMID:37138676; DOI:10.1016/j.jhepr.2023.100726"
)
t = selected_text(gse)
if not re.search(r"\bNC[123]\b", t, re.I):
    raise RuntimeError("GSE207281: NC control samples not detected.")
if not re.search(r"\bHFD[123]\b", t, re.I):
    raise RuntimeError("GSE207281: HFD control-genotype samples not detected.")
if re.search(r"\bNKO[123]\b|\bHKO[123]\b", t, re.I):
    raise RuntimeError("GSE207281: Pld1 knockout samples detected in frozen primary contrast.")

fill_var(idx, "sex", "MALE", "PUBLICATION_METHODS", citation,
         "Male Pld1-Flox and hepatocyte-specific Pld1-KO mice were used.")
fill_var(idx, "hfd_duration_weeks", "20", "PUBLICATION_METHODS", citation,
         "The HFD-induced NAFLD model was maintained for 20 weeks.")
fill_var(idx, "hfd_fat_kcal_percent", "60", "PUBLICATION_METHODS", citation,
         "HFD was Research Diets D12492, 60 kcal% fat.")
fill_var(idx, "control_diet_type", "NORMAL_DIET", "PUBLICATION_METHODS", citation,
         "Comparator was explicitly termed normal control diet (NCD).")
mark_not_reported(
    idx, "strain", "PUBLICATION_METHODS", citation,
    "The accessible Methods do not explicitly state the RNA-seq cohort C57BL/6 J/N substrain/background."
)
x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected samples are NC1-3 versus HFD1-3; NKO/HKO knockout samples excluded."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH05_PARTIAL_REVIEWED"
attempts.append({
    "GSE": "GSE207281",
    "unresolved_fields": "strain",
    "reason": "Publication does not explicitly state the RNA-seq cohort substrain/background."
})

# GSE282595 — PARTIAL
gse = "GSE282595"
idx = set_pub(
    gse, "40680014", "10.2337/db24-1173",
    "Concerted Actions of FoxO1 and PPARα in Hepatic Gene Expression and Metabolic Adaptation",
)
citation = (
    "Kitamoto et al. Diabetes. 2025;74:1748-1760. "
    "PMID:40680014; DOI:10.2337/db24-1173"
)
t = selected_text(gse)
if not re.search(r"RD_WT_Fast", t, re.I):
    raise RuntimeError("GSE282595: RD_WT_Fast controls not detected.")
if not re.search(r"HFD_WT_Fast", t, re.I):
    raise RuntimeError("GSE282595: HFD_WT_Fast samples not detected.")
if re.search(r"LFKO|LPKO|LFPKO", t, re.I):
    raise RuntimeError("GSE282595: liver-specific KO samples detected in frozen primary contrast.")

fill_var(idx, "sex", "MALE", "PUBLICATION_METHODS", citation,
         "Animals and diets section explicitly states male mice.")
fill_var(idx, "hfd_duration_weeks", "3", "GEO_GSM",
         "GSE282595 selected WT GSM treatment protocol + publication",
         "Selected WT HFD and regular-diet GSMs explicitly state diet was given for 3 weeks before tissue harvest.")
fill_var(idx, "hfd_fat_kcal_percent", "60", "PUBLICATION_METHODS", citation,
         "HFD was Research Diets D12492i, 60% kcal.")
fill_var(idx, "control_diet_type", "CHOW", "PUBLICATION_METHODS", citation,
         "Control was standard chow, PicoLab Rodent Diet 20 (5053; Purina Mills).")
mark_not_reported(
    idx, "strain", "PUBLICATION_METHODS", citation,
    "The paper identifies primary hepatocyte WT preparations as C57BL/6J, but does not explicitly map the selected WT floxed RNA-seq cohort to a J/N substrain."
)
x.loc[idx, "sample_role_identity_status"] = "MANUAL_CONFIRMED"
x.loc[idx, "sample_role_identity_evidence"] = (
    "Frozen selected samples are RD_WT_Fast vs HFD_WT_Fast; LFKO/LPKO/LFPKO samples excluded."
)
x.loc[idx, "overall_adjudication_status"] = "BATCH05_PARTIAL_REVIEWED"
x.loc[idx, "general_note"] = "Selected RNA-seq samples are 4-hour fasted."
attempts.append({
    "GSE": "GSE282595",
    "unresolved_fields": "strain",
    "reason": "No explicit J/N mapping for the selected WT floxed RNA-seq cohort."
})

# QA
batch_gses = ["GSE247985", "GSE137365", "GSE207281", "GSE282595"]
b = x[x["GSE"].isin(batch_gses)].copy()
assert b["GSE"].nunique() == 4
assert b["Contrast_id"].nunique() == 4

for gse in ["GSE247985", "GSE137365"]:
    z = b[b["GSE"].eq(gse)]
    for var in ["sex", "strain", "hfd_duration_weeks", "hfd_fat_kcal_percent", "control_diet_type"]:
        assert (z[f"{var}_final"] != "").all(), (gse, var)
    assert set(z["overall_adjudication_status"]) == {"BATCH05_ADJUDICATED"}

for gse in ["GSE207281", "GSE282595"]:
    z = b[b["GSE"].eq(gse)]
    assert set(z["strain_final"]) == {""}
    for var in ["sex", "hfd_duration_weeks", "hfd_fat_kcal_percent", "control_diet_type"]:
        assert (z[f"{var}_final"] != "").all(), (gse, var)
    assert set(z["overall_adjudication_status"]) == {"BATCH05_PARTIAL_REVIEWED"}

assert set(b[b["GSE"].eq("GSE247985")]["hfd_duration_weeks_final"]) == {"10"}
assert set(b[b["GSE"].eq("GSE137365")]["hfd_duration_weeks_final"]) == {"16"}
assert set(b[b["GSE"].eq("GSE137365")]["hfd_fat_kcal_percent_final"]) == {"39.5"}
assert set(b[b["GSE"].eq("GSE207281")]["hfd_duration_weeks_final"]) == {"20"}
assert set(b[b["GSE"].eq("GSE282595")]["hfd_duration_weeks_final"]) == {"3"}

# Preserve all previously reviewed statuses
prev = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
old_reviewed = prev[
    prev["overall_adjudication_status"].str.contains(
        "ADJUDICATED|PARTIAL_REVIEWED", regex=True, na=False
    )
]["GSE"].unique()

for gse in old_reviewed:
    a = prev[prev["GSE"].eq(gse)].sort_values("Contrast_id")
    bb = x[x["GSE"].eq(gse)].sort_values("Contrast_id")
    assert a["overall_adjudication_status"].tolist() == bb["overall_adjudication_status"].tolist(), gse

x.to_csv(OUT, sep="\t", index=False)
pd.DataFrame(attempts).to_csv(ATTEMPT_OUT, sep="\t", index=False)

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

print("===== HFD MODERATOR MANUAL ADJUDICATION BATCH05 PASS =====")
print("Reviewed GSE this batch: 4")
print("Fully adjudicated this batch: 2")
print("Partially reviewed this batch: 2")
print("Cumulative fully adjudicated GSE:", int((sdf["status"] == "FULLY_ADJUDICATED").sum()))
print("Cumulative partial-reviewed GSE:", int((sdf["status"] == "PARTIAL_REVIEWED").sum()))
print("Cumulative not-yet-reviewed GSE:", int((sdf["status"] == "NOT_YET_REVIEWED").sum()))

print("\n===== BATCH05 FINAL STATUS =====")
cols = [
    "GSE", "Contrast_id", "overall_adjudication_status",
    "sex_final", "strain_final", "hfd_duration_weeks_final",
    "hfd_fat_kcal_percent_final", "control_diet_type_final",
]
print(b[cols].sort_values(["GSE","Contrast_id"]).to_string(index=False))

print("\n===== SOURCE-INSUFFICIENT FIELDS =====")
print(pd.DataFrame(attempts).to_string(index=False))

print("\nKEY CORRECTIONS:")
print("- GSE137365: auto HFD fat=5% was D-allulose supplementation, not fat%; final=39.5%.")
print("- GSE247985: liver RNA-seq is the 10-week male cohort; 16-week extension is separate.")
print("- GSE207281 and GSE282595: duration upgraded to HIGH; strain remains unresolved.")

print("\nSaved:")
print(OUT.relative_to(ROOT))
print(ATTEMPT_OUT.relative_to(ROOT))

print("\nIMPORTANT:")
print("- No strain/substrain was inferred where source mapping was insufficient.")
print("- No intervention/genotype stratum was mixed into frozen primary contrasts.")
print("- No frozen eligibility or contrast membership changed.")
