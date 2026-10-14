#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd

# ============================================================
# HFD moderator manual adjudication — Batch10 v1
#
# Canonical input:
#   HFD_primary105_moderator_manual_adjudication_batch09_v1.tsv
#
# Targets:
#   GSE121343
#     - HFD fat = 54 energy% (source reports 54 kJ% fat; energy fraction
#       is unit-invariant between kJ and kcal)
#     - sex = MALE
#     - control diet = LFD
#     - duration intentionally NOT changed because public GEO/publication
#       timeline descriptions are not fully concordant.
#
#   GSE124463
#     - HFD fat = 60 kcal%
#     - duration = 20 wk
#
#   GSE176112
#     - HFD fat = 60 kcal%
#     - duration = 8 wk
#     - sex = MALE
#     - control diet = NORMAL_DIET
#
# Frozen biological/technical eligibility, GSM membership, contrast
# membership, sample-role identity, GEO-auto evidence, and untargeted
# moderator fields are immutable.
# ============================================================

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=__import__("os").environ.get("HFD_PROJECT_ROOT", "."),
)
args = ap.parse_args()

ROOT = Path(args.root)
OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"

MASTER = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch09_v1.tsv"
OUT = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch10_v1.tsv"
AUDIT_OUT = OUTDIR / "HFD_fat_manual_adjudication_batch10_v1_audit.tsv"
NOTE_OUT = OUTDIR / "HFD_batch10_v1_fat_adjudication_note.txt"

if not MASTER.exists():
    raise FileNotFoundError(MASTER)
if OUT.exists():
    raise FileExistsError(
        f"Refusing to overwrite existing Batch10 master: {OUT}\n"
        "Use batch10_v1_1 or later for corrections."
    )

x = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
prev = x.copy(deep=True)

assert x["GSE"].nunique() == 64
assert x["Contrast_id"].nunique() == 105

TARGETS = {
    "GSE121343": {"hfd_fat_kcal_percent", "sex", "control_diet_type"},
    "GSE124463": {"hfd_fat_kcal_percent", "hfd_duration_weeks"},
    "GSE176112": {"hfd_fat_kcal_percent", "hfd_duration_weeks", "sex", "control_diet_type"},
}

def one(gse):
    z = x[x["GSE"].eq(gse)]
    if len(z) != 1:
        raise RuntimeError(f"{gse}: expected one frozen contrast row; observed {len(z)}")
    return z.index[0]

idx = {gse: one(gse) for gse in TARGETS}

def effective(row, var):
    if row[f"{var}_adjudication_status"] != "NOT_REVIEWED":
        return (
            row[f"{var}_final"].strip(),
            row[f"{var}_final_confidence"].strip(),
            "MANUAL",
        )
    if row[f"{var}_geo_auto_status"] == "RESOLVED":
        return (
            row[f"{var}_geo_auto_value"].strip(),
            row[f"{var}_geo_auto_confidence"].strip(),
            "GEO_AUTO",
        )
    return "", "", "NONE"

def numeric_ready(val):
    try:
        z = float(val)
        return 0 < z <= 100
    except Exception:
        return False

def strict_fat(row):
    val, conf, origin = effective(row, "hfd_fat_kcal_percent")
    if not numeric_ready(val):
        return False
    if origin == "MANUAL":
        return conf == "HIGH"
    if origin == "GEO_AUTO":
        return conf in {"MEDIUM", "HIGH"}
    return False

# Canonical Batch09 anchor.
strict_fat_gses = {
    r["GSE"] for _, r in x.iterrows() if strict_fat(r)
}
assert len(strict_fat_gses) == 36, len(strict_fat_gses)

# Target fields must still be unresolved/pending before modification.
for gse, vars_ in TARGETS.items():
    i = idx[gse]
    for var in vars_:
        val, conf, origin = effective(x.loc[i], var)
        assert not val, (gse, var, val, conf, origin)
        assert x.loc[i, f"{var}_adjudication_status"] == "NOT_REVIEWED", (
            gse, var, x.loc[i, f"{var}_adjudication_status"]
        )

# GSE121343 duration is explicitly protected from modification.
i121 = idx["GSE121343"]
duration121_before = prev.loc[i121, [
    c for c in x.columns if c.startswith("hfd_duration_weeks_")
]].copy()

def fill_var(i, var, value, source_type, citation, evidence,
             confidence="HIGH", status="MANUAL_CONFIRMED", final_note=""):
    x.loc[i, f"{var}_manual_value"] = str(value)
    x.loc[i, f"{var}_manual_source_type"] = source_type
    x.loc[i, f"{var}_manual_source_citation"] = citation
    x.loc[i, f"{var}_manual_evidence"] = evidence
    x.loc[i, f"{var}_manual_confidence"] = confidence
    x.loc[i, f"{var}_adjudication_status"] = status
    x.loc[i, f"{var}_final"] = str(value)
    x.loc[i, f"{var}_final_confidence"] = confidence
    x.loc[i, f"{var}_final_source"] = source_type
    x.loc[i, f"{var}_final_note"] = final_note

def set_pub_if_safe(i, pmid, doi, title):
    for col, new in [
        ("publication_pmid", pmid),
        ("publication_doi", doi),
        ("publication_title", title),
    ]:
        old = x.loc[i, col]
        if old and new and old != new:
            raise RuntimeError(
                f"{x.loc[i,'GSE']}: conflicting {col}: existing={old!r}, new={new!r}"
            )
        if not old:
            x.loc[i, col] = new

# ============================================================
# GSE121343
# Pradas-Juni et al. Nat Commun. 2020.
# PMID 32005828; DOI 10.1038/s41467-020-14323-y
#
# Direct selected GSM treatment metadata:
# - HFD D12492(I) mod (Sniff)
# - 54 kJ% fat
# - CD D12450B* mod LS (Sniff), 11 kJ% fat
#
# Primary Methods identify the LNA cohort as male.
#
# Important: 54 kJ% is a percent of dietary energy. Percent energy is
# dimensionless and unchanged by expressing energy in kJ or kcal.
# This is NOT a wt% -> kcal% conversion.
#
# Duration intentionally left untouched because GEO and final paper
# descriptions differ in age/injection/diet timeline details.
# ============================================================
gse = "GSE121343"
i = idx[gse]
citation_geo = "GEO:GSE121343 selected CTRL/HFD LNA-cohort GSM treatment metadata"
citation_pub = (
    "Pradas-Juni et al. Nat Commun. 2020;11:644. "
    "PMID:32005828; DOI:10.1038/s41467-020-14323-y; GEO:GSE121343"
)
fill_var(
    i, "hfd_fat_kcal_percent", "54",
    "GEO_GSM+PUBLICATION_METHODS", citation_pub,
    "Selected GSE121343 GSM treatment metadata and the primary Methods identify "
    "the modified HFD as D12492(I) mod (Sniff) containing 54 kJ% fat. Because "
    "this is explicitly a percentage of dietary energy, the fraction is "
    "unit-invariant between kJ and kcal; final HFD fat energy percentage is 54. "
    "No ingredient weight-percent conversion is performed."
)
fill_var(
    i, "sex", "MALE",
    "PUBLICATION_METHODS", citation_pub,
    "The primary Methods explicitly describe the in vivo LNA cohort underlying "
    "GSE121343 as male C57BL/6N mice."
)
fill_var(
    i, "control_diet_type", "LFD",
    "GEO_GSM+PUBLICATION_METHODS", citation_pub,
    "The matched controlled-feeding comparator is D12450B* mod LS (Sniff), "
    "reported as 11% of dietary energy from fat. It is therefore a directly "
    "defined low-fat control diet (LFD)."
)
set_pub_if_safe(
    i, "32005828", "10.1038/s41467-020-14323-y",
    "A MAFG-lncRNA axis links systemic nutrient abundance to hepatic glucose metabolism",
)

# ============================================================
# GSE124463
# Qin et al. Nat Commun. 2020.
# PMID 32075973; DOI 10.1038/s41467-020-14796-x
#
# Primary article + GEO:
# - male C57BL/6
# - D12492 lipid-rich diet
# - 60% of caloric content from lipid
# - 20 weeks
# - liver RNA-seq deposited as GSE124463
# ============================================================
gse = "GSE124463"
i = idx[gse]
citation = (
    "Qin et al. Nat Commun. 2020;11:962. "
    "PMID:32075973; DOI:10.1038/s41467-020-14796-x; GEO:GSE124463"
)
fill_var(
    i, "hfd_fat_kcal_percent", "60",
    "PUBLICATION_METHODS", citation,
    "The primary article directly linked to GSE124463 states that the "
    "lipid-rich diet was Research Diets D12492 with 60% of caloric content "
    "from lipid."
)
fill_var(
    i, "hfd_duration_weeks", "20",
    "PUBLICATION_METHODS", citation,
    "Male C57BL/6 mice were placed on D12492 lipid-rich diet or the matched "
    "carbohydrate-rich diet for 20 weeks; the liver RNA-seq is deposited as "
    "GSE124463."
)
set_pub_if_safe(
    i, "32075973", "10.1038/s41467-020-14796-x",
    "Alterations in promoter interaction landscape and transcriptional network underlying metabolic adaptation to diet",
)

# ============================================================
# GSE176112
# Sun et al. EMBO Rep. 2021.
# PMID 34459531; DOI 10.15252/embr.202152537
#
# Primary Methods:
# - p53+/+ and p53-/- C57BL/6N male mice
# - normal diet or HFD D12492
# - HFD = 60% kcal from fat
# - 8 weeks
# - liver RNA-seq deposited as GSE176112
# ============================================================
gse = "GSE176112"
i = idx[gse]
citation = (
    "Sun et al. EMBO Rep. 2021;22:e52537. "
    "PMID:34459531; DOI:10.15252/embr.202152537; GEO:GSE176112"
)
fill_var(
    i, "hfd_fat_kcal_percent", "60",
    "PUBLICATION_METHODS", citation,
    "The primary Methods directly identify the HFD as Research Diets D12492, "
    "60% kcal from fat (5.24 kcal/g), in the liver RNA-seq cohort."
)
fill_var(
    i, "hfd_duration_weeks", "8",
    "PUBLICATION_METHODS", citation,
    "For the HFD-induced obese mouse model used for liver RNA-seq, male "
    "p53 wild-type and knockout mice were fed normal diet or D12492 HFD for "
    "8 weeks before liver analysis."
)
fill_var(
    i, "sex", "MALE",
    "PUBLICATION_METHODS", citation,
    "The primary Methods explicitly identify the normal/HFD p53 wild-type "
    "and knockout RNA-seq cohort as male mice."
)
fill_var(
    i, "control_diet_type", "NORMAL_DIET",
    "PUBLICATION_METHODS", citation,
    "The matched comparator for the RNA-seq experiment is explicitly described "
    "as normal diet. No more specific chow-versus-purified-LFD classification "
    "is inferred."
)
set_pub_if_safe(
    i, "34459531", "10.15252/embr.202152537",
    "p53 transcriptionally regulates SQLE to repress cholesterol synthesis and tumor growth",
)

# ------------------------------------------------------------
# Immutability QA
# ------------------------------------------------------------

# Protect GSE121343 duration byte-for-byte.
duration_cols = [c for c in x.columns if c.startswith("hfd_duration_weeks_")]
assert prev.loc[i121, duration_cols].equals(x.loc[i121, duration_cols])

# GEO-auto fields are immutable everywhere.
geo_auto_cols = [c for c in x.columns if "_geo_auto_" in c]
assert prev[geo_auto_cols].equals(x[geo_auto_cols])

# Frozen identity and membership fields immutable everywhere.
frozen_cols = [
    "GSE", "Contrast_id", "n_CONTROL_GSM", "n_HFD_GSM",
    "manual_priority_score", "manual_review_issues",
    "sample_role_identity_status",
    "sample_role_identity_evidence",
    "sample_role_identity_note",
]
for c in frozen_cols:
    if c in x.columns:
        assert prev[c].equals(x[c]), c

pub_cols = {"publication_pmid", "publication_doi", "publication_title"}
allowed_by_gse = {
    gse: {
        c
        for var in vars_
        for c in x.columns
        if c.startswith(var + "_") and "_geo_auto_" not in c
    } | pub_cols
    for gse, vars_ in TARGETS.items()
}

changed = []
for row_i in x.index:
    gse = x.loc[row_i, "GSE"]
    for col in x.columns:
        if prev.loc[row_i, col] != x.loc[row_i, col]:
            changed.append((gse, x.loc[row_i, "Contrast_id"], col))
            assert gse in TARGETS, ("Unexpected GSE modified", gse, col)
            assert col in allowed_by_gse[gse], ("Unexpected field modified", gse, col)

assert {g for g, _, _ in changed} == set(TARGETS)

# All untargeted GSE rows byte-identical.
untargeted = ~x["GSE"].isin(TARGETS)
assert prev.loc[untargeted].reset_index(drop=True).equals(
    x.loc[untargeted].reset_index(drop=True)
)

expected = {
    "GSE121343": {
        "hfd_fat_kcal_percent": ("54", "HIGH"),
        "sex": ("MALE", "HIGH"),
        "control_diet_type": ("LFD", "HIGH"),
    },
    "GSE124463": {
        "hfd_fat_kcal_percent": ("60", "HIGH"),
        "hfd_duration_weeks": ("20", "HIGH"),
    },
    "GSE176112": {
        "hfd_fat_kcal_percent": ("60", "HIGH"),
        "hfd_duration_weeks": ("8", "HIGH"),
        "sex": ("MALE", "HIGH"),
        "control_diet_type": ("NORMAL_DIET", "HIGH"),
    },
}
for gse, vals in expected.items():
    i = idx[gse]
    for var, (v, conf) in vals.items():
        assert x.loc[i, f"{var}_final"] == v, (gse, var)
        assert x.loc[i, f"{var}_final_confidence"] == conf, (gse, var)

x.to_csv(OUT, sep="\t", index=False)

audit_rows = []
for gse, vars_ in TARGETS.items():
    i = idx[gse]
    for var in sorted(vars_):
        audit_rows.append({
            "GSE": gse,
            "Contrast_id": x.loc[i, "Contrast_id"],
            "Moderator": var,
            "Prior_final_value": prev.loc[i, f"{var}_final"],
            "Prior_final_confidence": prev.loc[i, f"{var}_final_confidence"],
            "Batch10_final_value": x.loc[i, f"{var}_final"],
            "Batch10_final_confidence": x.loc[i, f"{var}_final_confidence"],
            "Source_type": x.loc[i, f"{var}_manual_source_type"],
            "Source_citation": x.loc[i, f"{var}_manual_source_citation"],
            "Evidence": x.loc[i, f"{var}_manual_evidence"],
            "Action": "MANUAL_HIGH_DIRECT_SOURCE",
        })

audit_rows.append({
    "GSE": "GSE121343",
    "Contrast_id": x.loc[i121, "Contrast_id"],
    "Moderator": "hfd_duration_weeks",
    "Prior_final_value": prev.loc[i121, "hfd_duration_weeks_final"],
    "Prior_final_confidence": prev.loc[i121, "hfd_duration_weeks_final_confidence"],
    "Batch10_final_value": x.loc[i121, "hfd_duration_weeks_final"],
    "Batch10_final_confidence": x.loc[i121, "hfd_duration_weeks_final_confidence"],
    "Source_type": "BATCH10_REVIEW_NO_CHANGE",
    "Source_citation": citation_pub + "; GEO:GSE121343 selected GSM metadata",
    "Evidence": (
        "Duration was intentionally not adjudicated in Batch10 because GEO "
        "sample metadata and the final publication describe the LNA/diet timeline "
        "with non-identical age and treatment-duration details. No start-to-collection "
        "duration is inferred from the conflicting timeline descriptions."
    ),
    "Action": "DEFERRED_TIMELINE_RECONCILIATION",
})

pd.DataFrame(audit_rows).to_csv(AUDIT_OUT, sep="\t", index=False)

NOTE_OUT.write_text(
f"""HFD moderator manual adjudication — Batch10 v1
Canonical input: {MASTER.name}
Output: {OUT.name}

Direct HIGH upgrades

GSE121343
- HFD fat: 54 energy% -> final hfd_fat_kcal_percent = 54 HIGH.
  Source explicitly reports 54 kJ% fat. This is a percentage of energy,
  so the numeric fraction is unchanged in kcal units. This is NOT a
  weight-percent conversion.
- sex: MALE HIGH.
- control diet: LFD HIGH (D12450B* mod LS, 11 energy% fat).
- duration: intentionally unchanged/deferred because GEO and publication
  timeline descriptions are not fully concordant.

GSE124463
- HFD fat: 60 kcal% HIGH.
- duration: 20 wk HIGH.
- Primary article directly links GSE124463 and D12492 60% caloric lipid diet.

GSE176112
- HFD fat: 60 kcal% HIGH.
- duration: 8 wk HIGH.
- sex: MALE HIGH.
- control diet: NORMAL_DIET HIGH.
- Primary publication directly links liver RNA-seq to GSE176112.

Frozen/immutable
- biological eligibility unchanged
- technical eligibility unchanged
- GSM membership unchanged
- contrast membership unchanged
- sample-role identity unchanged
- all GEO-auto fields unchanged
- all untargeted GSE rows byte-identical

Expected independent-GSE strict coverage after Batch10
- HFD duration: 40
- sex: 42
- strain: 41
- HFD fat kcal%: 39
- control diet: 52

Expected independent-GSE inclusive coverage
- HFD duration: 42
- HFD fat kcal%: 40
- control diet: 54

Versioning
- Batch09 v1 is not overwritten.
- Batch10 v1 refuses overwrite.
- Corrections must supersede as Batch10 v1.1 or later.
""",
    encoding="utf-8",
)

print("===== HFD MODERATOR MANUAL ADJUDICATION BATCH10 v1 =====")
print("Canonical input:", MASTER.relative_to(ROOT))
print("Target GSE:", ", ".join(TARGETS))
print("Changed cells:", len(changed))
print("\nFinal target values:")
show = []
for gse in TARGETS:
    i = idx[gse]
    show.append({
        "GSE": gse,
        "sex": x.loc[i, "sex_final"],
        "strain": x.loc[i, "strain_final"],
        "duration": x.loc[i, "hfd_duration_weeks_final"],
        "fat_kcal_pct": x.loc[i, "hfd_fat_kcal_percent_final"],
        "control": x.loc[i, "control_diet_type_final"],
    })
print(pd.DataFrame(show).to_string(index=False))
print("\nGSE121343 duration unchanged:",
      repr(x.loc[i121, "hfd_duration_weeks_final"]),
      repr(x.loc[i121, "hfd_duration_weeks_final_confidence"]))
print("\nSaved:")
print(" ", OUT.relative_to(ROOT))
print(" ", AUDIT_OUT.relative_to(ROOT))
print(" ", NOTE_OUT.relative_to(ROOT))
print("QA PASS: frozen fields and all untargeted GSE rows are unchanged.")
