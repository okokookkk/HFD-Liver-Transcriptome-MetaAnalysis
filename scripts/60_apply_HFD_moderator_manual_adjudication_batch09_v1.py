#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd

# ============================================================
# HFD moderator manual adjudication — Batch09 v1
#
# Canonical input:
#   HFD_primary105_moderator_manual_adjudication_batch08_v1.tsv
#
# Primary goal:
#   Upgrade directly source-supported HFD fat kcal% values.
#
# Batch09 target GSE:
#   GSE124777
#   GSE137022
#   GSE145665
#   GSE146627
#   GSE193785
#
# Opportunistic co-adjudication is limited to moderators directly
# supported by the same primary/GEO evidence:
#   GSE124777: fat%, duration, sex
#   GSE137022: fat%, duration, strain
#   GSE145665: fat%, duration, sex
#   GSE146627: fat%, duration
#   GSE193785: fat%, duration, sex
#
# GSE145840 is intentionally NOT changed:
#   existing 60 wk? NO — existing fat=60 MEDIUM remains non-strict
#   because a direct primary/GEO kcal% anchor sufficient for HIGH was
#   not established in this Batch09 pass.
#
# Frozen biological/technical eligibility, GSM membership, contrast
# membership, sample-role identity, all GEO-auto fields, and all
# untargeted moderator fields are immutable.
# ============================================================

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=__import__("os").environ.get("HFD_PROJECT_ROOT", "."),
)
args = ap.parse_args()

ROOT = Path(args.root)
OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"

MASTER = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch08_v1.tsv"
OUT = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch09_v1.tsv"
AUDIT_OUT = OUTDIR / "HFD_fat_manual_adjudication_batch09_v1_audit.tsv"
NOTE_OUT = OUTDIR / "HFD_batch09_v1_fat_adjudication_note.txt"

if not MASTER.exists():
    raise FileNotFoundError(MASTER)
if OUT.exists():
    raise FileExistsError(
        f"Refusing to overwrite existing Batch09 master: {OUT}\n"
        "Use batch09_v1_1 or later for any correction."
    )

x = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
prev = x.copy(deep=True)

assert x["GSE"].nunique() == 64
assert x["Contrast_id"].nunique() == 105

TARGETS = {
    "GSE124777": {"hfd_fat_kcal_percent", "hfd_duration_weeks", "sex"},
    "GSE137022": {"hfd_fat_kcal_percent", "hfd_duration_weeks", "strain"},
    "GSE145665": {"hfd_fat_kcal_percent", "hfd_duration_weeks", "sex"},
    "GSE146627": {"hfd_fat_kcal_percent", "hfd_duration_weeks"},
    "GSE193785": {"hfd_fat_kcal_percent", "hfd_duration_weeks", "sex"},
}

def one(gse):
    z = x[x["GSE"].eq(gse)]
    if len(z) != 1:
        raise RuntimeError(f"{gse}: expected exactly one frozen contrast row, observed {len(z)}")
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

# Canonical Batch08 anchor.
strict_fat_gses = {
    r["GSE"] for _, r in x.iterrows() if strict_fat(r)
}
assert len(strict_fat_gses) == 31, len(strict_fat_gses)

# All five Batch09 fat targets came from PENDING_MISSING in the
# canonical Batch08 queue: they must not already have an effective fat%.
for gse, i in idx.items():
    val, conf, origin = effective(x.loc[i], "hfd_fat_kcal_percent")
    assert not numeric_ready(val), (gse, val, conf, origin)
    assert x.loc[i, "hfd_fat_kcal_percent_adjudication_status"] == "NOT_REVIEWED", (
        gse, x.loc[i, "hfd_fat_kcal_percent_adjudication_status"]
    )

# Expected opportunistic fields must also be unresolved/pending,
# preventing silent overwrite of a newer server-side adjudication.
for gse, var in [
    ("GSE124777", "hfd_duration_weeks"),
    ("GSE124777", "sex"),
    ("GSE137022", "hfd_duration_weeks"),
    ("GSE137022", "strain"),
    ("GSE145665", "hfd_duration_weeks"),
    ("GSE145665", "sex"),
    ("GSE146627", "hfd_duration_weeks"),
    ("GSE193785", "hfd_duration_weeks"),
    ("GSE193785", "sex"),
]:
    i = idx[gse]
    val, conf, origin = effective(x.loc[i], var)
    assert not val, (gse, var, val, conf, origin)
    assert x.loc[i, f"{var}_adjudication_status"] == "NOT_REVIEWED", (
        gse, var, x.loc[i, f"{var}_adjudication_status"]
    )

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
    # Do not silently overwrite conflicting existing publication metadata.
    for col, new in [
        ("publication_pmid", pmid),
        ("publication_doi", doi),
        ("publication_title", title),
    ]:
        old = x.loc[i, col]
        if old and new and old != new:
            raise RuntimeError(
                f"{x.loc[i, 'GSE']}: conflicting {col}: existing={old!r}, new={new!r}"
            )
        if not old:
            x.loc[i, col] = new

# ============================================================
# GSE124777
# Shon et al., Nutrients 2020.
# PMID 33003339; DOI 10.3390/nu12102973
#
# Primary publication:
# - male C57BL/6J
# - HFD 45 kcal% fat
# - 16 weeks
# - liver RNA-seq deposited as GSE124777
# Frozen contrast membership is NOT re-screened/changed.
# ============================================================
gse = "GSE124777"
i = idx[gse]
citation = (
    "Shon et al. Nutrients. 2020;12(10):2973. "
    "PMID:33003339; DOI:10.3390/nu12102973; GEO:GSE124777"
)
fill_var(
    i, "hfd_fat_kcal_percent", "45",
    "PUBLICATION_METHODS", citation,
    "The primary publication directly links its liver RNA-seq to GSE124777 and "
    "states that the HFD group received a diet containing 45% kcal from fat."
)
fill_var(
    i, "hfd_duration_weeks", "16",
    "PUBLICATION_METHODS", citation,
    "Male C57BL/6J mice were fed normal diet, HFD, HFD+CLE, or HFD+LU for "
    "16 weeks; at the end of the experimental period liver was collected and "
    "the liver RNA-seq data were deposited as GSE124777."
)
fill_var(
    i, "sex", "MALE",
    "PUBLICATION_METHODS", citation,
    "The primary publication explicitly states that the GSE124777 mouse cohort "
    "used male C57BL/6J mice."
)
set_pub_if_safe(
    i, "33003339", "10.3390/nu12102973",
    "Plasma Lipidomics Reveals Insights into Anti-Obesity Effect of Chrysanthemum morifolium Ramat Leaves and Its Constituent Luteolin in High-Fat Diet-Induced Dyslipidemic Mice",
)

# ============================================================
# GSE137022
# FAM13A affects body fat distribution and adipocyte function.
# Nature Communications 2020.
# PMID 32193374; DOI 10.1038/s41467-020-15291-z
#
# Primary publication:
# - WT/KO developed on C57BL6 background
# - HFD Research Diets D12451 = Fat 45% kcal
# - HFD from age 8 wk for 14 wk, killed at age 22 wk
# - RNA-seq deposited as GSE137022
# ============================================================
gse = "GSE137022"
i = idx[gse]
citation = (
    "FAM13A affects body fat distribution and adipocyte function. "
    "Nat Commun. 2020;11:1465. PMID:32193374; "
    "DOI:10.1038/s41467-020-15291-z; GEO:GSE137022"
)
fill_var(
    i, "hfd_fat_kcal_percent", "45",
    "PUBLICATION_METHODS", citation,
    "The primary publication identifies the HFD as Research Diets D12451 "
    "and explicitly reports fat as 45% kcal."
)
fill_var(
    i, "hfd_duration_weeks", "14",
    "PUBLICATION_METHODS", citation,
    "Mice began HFD at 8 weeks of age and were fed HFD for 14 weeks until "
    "they were killed at 22 weeks of age; RNA-seq was deposited as GSE137022."
)
fill_var(
    i, "strain", "C57BL/6",
    "PUBLICATION_METHODS", citation,
    "The primary publication states that control and Fam13a KO mice were "
    "developed on a C57BL6 background. No J/N substrain is inferred."
)
set_pub_if_safe(
    i, "32193374", "10.1038/s41467-020-15291-z",
    "FAM13A affects body fat distribution and adipocyte function",
)

# ============================================================
# GSE145665
#
# Direct untreated HFD GSM:
#   GSM4323503 (liver, M1)
#   - high fat diet
#   - D12492, Research Diets
#   - 24 weeks
#   - GEO text reports 60.9% fat but does not label that number as kcal%.
#
# For the variable hfd_fat_kcal_percent we therefore use the exact
# product code rather than treating the unlabeled 60.9% as kcal%.
# Research Diets' D12492 formula = 60 kcal% fat.
#
# Follow-up primary paper explicitly links RNA-seq to GSE145665 and
# states the cohort consisted of male mice.
# PMID 36249786; DOI 10.3389/fphar.2022.946360
# ============================================================
gse = "GSE145665"
i = idx[gse]
citation_geo = "GEO:GSE145665; selected untreated HFD GSM4323503"
citation_pub = (
    "Kaempferol attenuates nonalcoholic steatohepatitis by regulating serum "
    "and liver bile acid metabolism. Front Pharmacol. 2022;13:946360. "
    "PMID:36249786; DOI:10.3389/fphar.2022.946360; GEO:GSE145665"
)
fill_var(
    i, "hfd_fat_kcal_percent", "60",
    "GEO_GSM+DIET_MANUFACTURER_CODE",
    citation_geo + "; Research Diets D12492 formula",
    "Untreated HFD GSM4323503 directly identifies the diet as D12492 "
    "(Research Diets) and reports 24 weeks of feeding. The GSM's 60.9% fat "
    "wording is not explicitly labeled as kcal%, so it is not copied as the "
    "kcal% moderator. The exact D12492 manufacturer formula specifies "
    "60 kcal% fat; final HFD fat kcal% is therefore 60."
)
fill_var(
    i, "hfd_duration_weeks", "24",
    "GEO_GSM", citation_geo,
    "Untreated frozen-group HFD GSM4323503 directly states that mice were "
    "fed the high-fat diet for 24 weeks."
)
fill_var(
    i, "sex", "MALE",
    "PUBLICATION_METHODS", citation_pub,
    "The publication directly linked to GSE145665 describes the mouse cohort "
    "as specific-pathogen-free male mice and identifies the RNA-seq dataset "
    "as GSE145665."
)
set_pub_if_safe(
    i, "36249786", "10.3389/fphar.2022.946360",
    "Kaempferol attenuates nonalcoholic steatohepatitis by regulating serum and liver bile acid metabolism",
)

# ============================================================
# GSE146627
# Luo et al. RNA Biology 2021.
# PMID 32967529; DOI 10.1080/15476286.2020.1827193
#
# Primary Methods:
# - male C57BL/6
# - HFD 45% kcal from fat; MD12032
# - 8 weeks feeding
# - sacrifice followed by liver isolation
# GEO links PMID 32967529 to GSE146627.
# ============================================================
gse = "GSE146627"
i = idx[gse]
citation = (
    "Luo et al. RNA Biol. 2021;18(6):863-874. PMID:32967529; "
    "DOI:10.1080/15476286.2020.1827193; GEO:GSE146627"
)
fill_var(
    i, "hfd_fat_kcal_percent", "45",
    "PUBLICATION_METHODS", citation,
    "The primary publication reports the HFD as MD12032 "
    "(Jiangsu Medicience) with 45% kcal from fat."
)
fill_var(
    i, "hfd_duration_weeks", "8",
    "PUBLICATION_METHODS", citation,
    "Eight-week-old male mice were fed control diet or HFD for 8 weeks and "
    "were then sacrificed for liver isolation; GEO GSE146627 links this "
    "publication to the deposited liver RNA-seq."
)
set_pub_if_safe(
    i, "32967529", "10.1080/15476286.2020.1827193",
    "Comprehensive analysis of the translatome reveals the relationship between the translational and transcriptional control in high fat diet-induced liver steatosis",
)

# ============================================================
# GSE193785
# Jin et al. Nutrients 2022.
# PMID 35215514; DOI 10.3390/nu14040865
#
# Primary Methods:
# - 4-week-old male C57BL/6J, 1-week acclimation
# - HFD D12492 for 5 weeks to induce obesity
# - HFD-induced animals then remain in HFD-based groups for a further
#   5-week treatment period; tissue collected at experiment end
# - total HFD start -> collection for untreated HFD comparator = 10 wk
# - D12492; manufacturer definition = 60 kcal% fat
# - data availability = GSE193785
# ============================================================
gse = "GSE193785"
i = idx[gse]
citation = (
    "Jin et al. Nutrients. 2022;14(4):865. PMID:35215514; "
    "DOI:10.3390/nu14040865; GEO:GSE193785"
)
fill_var(
    i, "hfd_fat_kcal_percent", "60",
    "PUBLICATION_METHODS+DIET_MANUFACTURER_CODE", citation + "; Research Diets D12492",
    "The primary publication identifies the HFD as Research Diets D12492. "
    "The exact D12492 formula specifies 60 kcal% fat."
)
fill_var(
    i, "hfd_duration_weeks", "10",
    "PUBLICATION_METHODS", citation,
    "Mice received D12492 HFD for 5 weeks to induce obesity, then the "
    "HFD-induced animals were assigned to HFD-based groups for a further "
    "5-week experimental period before tissue collection. Under the "
    "prespecified HFD-start-to-collection definition, the untreated HFD "
    "comparator therefore has 10 weeks total HFD exposure."
)
fill_var(
    i, "sex", "MALE",
    "PUBLICATION_METHODS", citation,
    "The primary publication explicitly states that four-week-old male "
    "C57BL/6J mice were used; the generated dataset is GSE193785."
)
set_pub_if_safe(
    i, "35215514", "10.3390/nu14040865",
    "Anti-Obesity Effects of Polymethoxyflavone-Rich Fraction from Jinkyool (Citrus sunki Hort. ex Tanaka) Leaf on Obese Mice Induced by High-Fat Diet",
)

# ------------------------------------------------------------
# GSE145840 — explicit no-change QA
# ------------------------------------------------------------
z145840_before = prev[prev["GSE"].eq("GSE145840")].copy()
assert len(z145840_before) == 1
r145840 = z145840_before.iloc[0]
assert float(r145840["hfd_fat_kcal_percent_final"]) == 60.0
assert r145840["hfd_fat_kcal_percent_final_confidence"] == "MEDIUM"

# ------------------------------------------------------------
# Immutability QA
# ------------------------------------------------------------

# GEO-auto fields must remain identical everywhere.
geo_auto_cols = [c for c in x.columns if "_geo_auto_" in c]
assert prev[geo_auto_cols].equals(x[geo_auto_cols])

# Frozen identity/contrast columns must remain identical everywhere.
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

# Allow only targeted moderator fields plus publication metadata.
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

assert changed
assert {g for g, _, _ in changed} == set(TARGETS)

# GSE145840 must remain byte-identical.
assert prev[prev["GSE"].eq("GSE145840")].reset_index(drop=True).equals(
    x[x["GSE"].eq("GSE145840")].reset_index(drop=True)
)

# Untargeted GSE rows must be byte-identical in full.
untargeted = ~x["GSE"].isin(TARGETS)
assert prev.loc[untargeted].reset_index(drop=True).equals(
    x.loc[untargeted].reset_index(drop=True)
)

# Target final-value QA.
expected = {
    "GSE124777": {
        "hfd_fat_kcal_percent": ("45", "HIGH"),
        "hfd_duration_weeks": ("16", "HIGH"),
        "sex": ("MALE", "HIGH"),
    },
    "GSE137022": {
        "hfd_fat_kcal_percent": ("45", "HIGH"),
        "hfd_duration_weeks": ("14", "HIGH"),
        "strain": ("C57BL/6", "HIGH"),
    },
    "GSE145665": {
        "hfd_fat_kcal_percent": ("60", "HIGH"),
        "hfd_duration_weeks": ("24", "HIGH"),
        "sex": ("MALE", "HIGH"),
    },
    "GSE146627": {
        "hfd_fat_kcal_percent": ("45", "HIGH"),
        "hfd_duration_weeks": ("8", "HIGH"),
    },
    "GSE193785": {
        "hfd_fat_kcal_percent": ("60", "HIGH"),
        "hfd_duration_weeks": ("10", "HIGH"),
        "sex": ("MALE", "HIGH"),
    },
}
for gse, vals in expected.items():
    i = idx[gse]
    for var, (v, conf) in vals.items():
        assert x.loc[i, f"{var}_final"] == v, (gse, var)
        assert x.loc[i, f"{var}_final_confidence"] == conf, (gse, var)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
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
            "Batch09_final_value": x.loc[i, f"{var}_final"],
            "Batch09_final_confidence": x.loc[i, f"{var}_final_confidence"],
            "Source_type": x.loc[i, f"{var}_manual_source_type"],
            "Source_citation": x.loc[i, f"{var}_manual_source_citation"],
            "Evidence": x.loc[i, f"{var}_manual_evidence"],
            "Action": "MANUAL_HIGH_DIRECT_SOURCE",
        })

# Record GSE145840 review without changing master.
audit_rows.append({
    "GSE": "GSE145840",
    "Contrast_id": r145840["Contrast_id"],
    "Moderator": "hfd_fat_kcal_percent",
    "Prior_final_value": r145840["hfd_fat_kcal_percent_final"],
    "Prior_final_confidence": r145840["hfd_fat_kcal_percent_final_confidence"],
    "Batch09_final_value": r145840["hfd_fat_kcal_percent_final"],
    "Batch09_final_confidence": r145840["hfd_fat_kcal_percent_final_confidence"],
    "Source_type": "BATCH09_REVIEW_NO_CHANGE",
    "Source_citation": "Primary publication/GEO reviewed; no new direct exact kcal% anchor accepted for HIGH.",
    "Evidence": (
        "The biological HFD cohort linkage is supported, but this Batch09 pass "
        "did not establish a primary/GEO exact kcal% source sufficient to "
        "supersede the existing 60 MEDIUM adjudication."
    ),
    "Action": "MEDIUM_RETAINED",
})

pd.DataFrame(audit_rows).to_csv(AUDIT_OUT, sep="\t", index=False)

NOTE_OUT.write_text(
f"""HFD moderator manual adjudication — Batch09 v1
Canonical input: {MASTER.name}
Output: {OUT.name}

Primary target: HFD fat kcal% strict upgrades.

Direct HIGH upgrades
- GSE124777: fat 45 kcal%; duration 16 wk; sex MALE.
- GSE137022: fat 45 kcal%; duration 14 wk; strain C57BL/6.
- GSE145665: fat 60 kcal%; duration 24 wk; sex MALE.
  Note: untreated HFD GSM directly identifies D12492 and 24 weeks.
  The GSM's 60.9% fat text is not itself labeled kcal%; D12492's
  manufacturer caloric composition (60 kcal% fat) is used for the
  hfd_fat_kcal_percent variable rather than copying 60.9.
- GSE146627: fat 45 kcal%; duration 8 wk.
- GSE193785: fat 60 kcal%; duration 10 wk; sex MALE.
  Duration is 5 wk HFD induction + a further 5 wk under HFD before
  final tissue collection; no pre-intervention-only duration is used.

Reviewed, no change
- GSE145840: existing fat 60 MEDIUM retained. No HIGH upgrade in Batch09.

Frozen/immutable
- biological eligibility unchanged
- technical eligibility unchanged
- GSM membership unchanged
- contrast membership unchanged
- sample-role identity unchanged
- all GEO-auto fields unchanged
- all untargeted GSE rows byte-identical

Expected independent-GSE strict coverage after recomputation
- HFD duration: 38
- sex: 40
- strain: 41
- HFD fat kcal%: 36
- control diet: 50

Expected duration inclusive coverage: 40 GSE.
Existing duration strict-gap queue (GSE285614, GSE324134) remains conceptually separate.

Versioning
- Batch08 v1 is not overwritten.
- Batch09 v1 refuses overwrite.
- Corrections must supersede as Batch09 v1.1 or later.
""",
    encoding="utf-8",
)

print("===== HFD MODERATOR MANUAL ADJUDICATION BATCH09 v1 =====")
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
print("\nGSE145840 unchanged:",
      r145840["hfd_fat_kcal_percent_final"],
      r145840["hfd_fat_kcal_percent_final_confidence"])
print("\nSaved:")
print(" ", OUT.relative_to(ROOT))
print(" ", AUDIT_OUT.relative_to(ROOT))
print(" ", NOTE_OUT.relative_to(ROOT))
print("QA PASS: frozen fields and all untargeted GSE rows are unchanged.")
