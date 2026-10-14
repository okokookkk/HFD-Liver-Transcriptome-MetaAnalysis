#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd

# ============================================================
# HFD moderator manual adjudication — Batch08 v1
#
# Canonical input:
#   HFD_primary105_moderator_manual_adjudication_batch07_v1.tsv
#
# Scope:
#   - Correct GSE285380 duration only:
#       prior 13.0446 wk / MEDIUM -> reviewed UNRESOLVED
#   - Re-review GSE285614 but DO NOT modify it:
#       13.0446 wk / MEDIUM remains unchanged
#   - GSE324134 MUST remain untouched:
#       NOT_REVIEWED + GEO_AUTO_LOW 12 wk
#
# Frozen biological/technical eligibility, GSM membership,
# contrast membership, sample-role identity, and all non-duration
# moderators are immutable in this script.
# ============================================================

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=__import__("os").environ.get("HFD_PROJECT_ROOT", "."),
)
args = ap.parse_args()

ROOT = Path(args.root)
OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"

MASTER = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch07_v1.tsv"
OUT = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch08_v1.tsv"
AUDIT_OUT = OUTDIR / "HFD_duration_strict_gap_adjudication_batch08_v1_audit.tsv"
NOTE_OUT = OUTDIR / "HFD_batch08_v1_duration_strict_gap_note.txt"

if not MASTER.exists():
    raise FileNotFoundError(MASTER)
if OUT.exists():
    raise FileExistsError(
        f"Refusing to overwrite existing Batch08 master: {OUT}\n"
        "Use a superseding version (e.g. batch08_v1_1) if correction is required."
    )

x = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
prev = x.copy(deep=True)

# ------------------------------------------------------------
# Global canonical-input QA
# ------------------------------------------------------------
assert x["GSE"].nunique() == 64, x["GSE"].nunique()
assert x["Contrast_id"].nunique() == 105, x["Contrast_id"].nunique()

required_cols = [
    "GSE", "Contrast_id",
    "overall_adjudication_status",
    "sample_role_identity_status",
    "hfd_duration_weeks_geo_auto_value",
    "hfd_duration_weeks_geo_auto_status",
    "hfd_duration_weeks_geo_auto_confidence",
    "hfd_duration_weeks_manual_value",
    "hfd_duration_weeks_manual_source_type",
    "hfd_duration_weeks_manual_source_citation",
    "hfd_duration_weeks_manual_evidence",
    "hfd_duration_weeks_manual_confidence",
    "hfd_duration_weeks_adjudication_status",
    "hfd_duration_weeks_final",
    "hfd_duration_weeks_final_confidence",
    "hfd_duration_weeks_final_source",
    "hfd_duration_weeks_final_note",
]
missing = [c for c in required_cols if c not in x.columns]
if missing:
    raise RuntimeError(f"Batch07 schema missing expected columns: {missing}")

def one(gse):
    z = x[x["GSE"].eq(gse)]
    if len(z) != 1:
        raise RuntimeError(f"{gse}: expected exactly one contrast row, observed {len(z)}")
    return z.index[0]

def assert_final(gse, var, value, confidence):
    z = x[x["GSE"].eq(gse)]
    vals = set(z[f"{var}_final"])
    confs = set(z[f"{var}_final_confidence"])
    assert vals == {value}, (gse, var, vals, value)
    assert confs == {confidence}, (gse, var, confs, confidence)

# Canonical Batch07 anchors.
# These deliberately distinguish the server canonical Batch07 from
# alternate/provisional Batch07 files created elsewhere.
assert_final("GSE197649", "hfd_duration_weeks", "10", "HIGH")
assert_final("GSE267123", "hfd_duration_weeks", "20", "HIGH")
assert_final("GSE283060", "hfd_duration_weeks", "16", "HIGH")
assert_final("GSE288077", "hfd_duration_weeks", "16", "HIGH")

z274 = x[x["GSE"].eq("GSE274914")]
assert z274["Contrast_id"].nunique() == 4
assert set(z274["hfd_duration_weeks_final"]) == {"10"}
assert set(z274["hfd_duration_weeks_final_confidence"]) == {"HIGH"}

z283 = x[x["GSE"].eq("GSE283060")]
assert set(z283["strain_final"]) == {"C57BL/6J"}
assert set(z283["strain_final_confidence"]) == {"HIGH"}
assert set(z283["control_diet_type_final"]) == {"CHOW"}
assert set(z283["control_diet_type_final_confidence"]) == {"HIGH"}

z324 = x[x["GSE"].eq("GSE324134")]
assert len(z324) == 1
r324 = z324.iloc[0]
assert r324["overall_adjudication_status"] == "NOT_REVIEWED"
assert r324["hfd_duration_weeks_final"] == ""
assert r324["hfd_duration_weeks_final_confidence"] == ""
assert r324["hfd_duration_weeks_adjudication_status"] == "NOT_REVIEWED"
assert r324["hfd_duration_weeks_geo_auto_status"] == "RESOLVED"
assert r324["hfd_duration_weeks_geo_auto_confidence"] == "LOW"
assert float(r324["hfd_duration_weeks_geo_auto_value"]) == 12.0

# Remaining strict-gap candidates must enter Batch08 with the expected
# Batch07 manual MEDIUM values.
i380 = one("GSE285380")
i614 = one("GSE285614")

for gse, i in [("GSE285380", i380), ("GSE285614", i614)]:
    assert float(x.loc[i, "hfd_duration_weeks_final"]) == 13.0446, (
        gse, x.loc[i, "hfd_duration_weeks_final"]
    )
    assert x.loc[i, "hfd_duration_weeks_final_confidence"] == "MEDIUM", (
        gse, x.loc[i, "hfd_duration_weeks_final_confidence"]
    )
    assert x.loc[i, "hfd_duration_weeks_adjudication_status"] != "NOT_REVIEWED", gse

# ------------------------------------------------------------
# GSE285380 — correction
#
# Evidence logic:
# - selected liver comparator is lean/no-estrogen vs obese/no-estrogen
# - GEO GSM treatment protocol gives chow/HFHSD for "3 months"
# - publication supplement separately describes a 5-week experimental
#   period followed by sacrifice/liver collection
# - source does not explicitly establish total HFD start->collection
#   exposure for the frozen no-estrogen obese comparator
#
# Therefore:
#   - withdraw 13.0446-wk MEDIUM as total HFD duration
#   - do NOT infer 3 months + 5 weeks
#   - final = NOT_REPORTED/UNRESOLVED
# ------------------------------------------------------------

citation_380 = (
    "Mohammad et al. J Clin Invest. 2025;135(22):e193976. "
    "PMID:41243971; DOI:10.1172/JCI193976; GEO:GSE285380; "
    "JCI Supplemental Methods"
)
evidence_380 = (
    "GEO GSM treatment metadata report normal chow or HFHSD feeding for 3 months. "
    "The JCI Supplemental Methods separately describe the four-group mouse experiment, "
    "a 5-week experimental period with weekly plasma sampling, and sacrifice with liver "
    "collection after the fifth week. The public sources do not explicitly establish "
    "that the stated 3 months equals the full HFD-start-to-liver-collection interval "
    "for the frozen no-estrogen obese comparator, nor do they explicitly state the "
    "total diet-exposure duration through the 5-week follow-up. Therefore neither "
    "13.0446 weeks nor 3 months + 5 weeks is accepted as total HFD duration."
)
note_380 = (
    "Batch08 correction: prior 13.0446-week MEDIUM value is superseded. "
    "The 3-month diet statement cannot be safely equated with total HFD "
    "start-to-collection exposure because liver collection follows a separately "
    "described 5-week experimental period. Exact total duration remains "
    "NOT_REPORTED/UNRESOLVED."
)

x.loc[i380, "hfd_duration_weeks_manual_value"] = ""
x.loc[i380, "hfd_duration_weeks_manual_source_type"] = "PUBLICATION_SUPPLEMENT+GEO_GSM"
x.loc[i380, "hfd_duration_weeks_manual_source_citation"] = citation_380
x.loc[i380, "hfd_duration_weeks_manual_evidence"] = evidence_380
x.loc[i380, "hfd_duration_weeks_manual_confidence"] = "UNRESOLVED"
x.loc[i380, "hfd_duration_weeks_adjudication_status"] = "NOT_REPORTED"
x.loc[i380, "hfd_duration_weeks_final"] = ""
x.loc[i380, "hfd_duration_weeks_final_confidence"] = ""
x.loc[i380, "hfd_duration_weeks_final_source"] = ""
x.loc[i380, "hfd_duration_weeks_final_note"] = note_380

# ------------------------------------------------------------
# GSE285614 — reviewed, no master modification
#
# Publication directly anchors the bulk-liver RNA-seq cohort to 3 months
# of chow/HFD exposure. Existing standardization = 13.0446 wk MEDIUM.
# No exact week/day source -> no HIGH upgrade.
# ------------------------------------------------------------

citation_614 = (
    "Guan et al. J Clin Invest. 2025;135(21):e190635. "
    "PMID:41178716; DOI:10.1172/JCI190635; GEO:GSE285614"
)
evidence_614 = (
    "Methods state that mice aged 8 weeks were placed on HFD (D12492, 60% kcal fat) "
    "or chow. The bulk liver RNA-seq experiment states that male Bcl2AdTG mice and "
    "WT littermates were fed chow or HFD for 3 months and liver tissues were subjected "
    "to RNA-seq. The biological HFD-start-to-liver-collection interval is therefore "
    "directly linked to the frozen RNA-seq cohort as 3 months, but no exact week/day "
    "duration is reported. Existing 13.0446-week MEDIUM standardization is retained."
)

# ------------------------------------------------------------
# Immutability QA
# ------------------------------------------------------------

# Only these 10 duration adjudication/final fields may change,
# and only on GSE285380__C01.
allowed_changed_cols = {
    "hfd_duration_weeks_manual_value",
    "hfd_duration_weeks_manual_source_type",
    "hfd_duration_weeks_manual_source_citation",
    "hfd_duration_weeks_manual_evidence",
    "hfd_duration_weeks_manual_confidence",
    "hfd_duration_weeks_adjudication_status",
    "hfd_duration_weeks_final",
    "hfd_duration_weeks_final_confidence",
    "hfd_duration_weeks_final_source",
    "hfd_duration_weeks_final_note",
}

changed = []
for idx in x.index:
    for col in x.columns:
        if prev.loc[idx, col] != x.loc[idx, col]:
            changed.append((idx, x.loc[idx, "GSE"], x.loc[idx, "Contrast_id"], col))
            assert idx == i380, (
                "Unexpected row modified", idx, x.loc[idx, "GSE"],
                x.loc[idx, "Contrast_id"], col
            )
            assert col in allowed_changed_cols, (
                "Unexpected column modified", x.loc[idx, "Contrast_id"], col
            )

assert changed, "No changes detected; expected GSE285380 duration correction."
assert {z[2] for z in changed} == {"GSE285380__C01"}

# Explicitly guarantee all non-duration fields are byte-equivalent.
non_duration_cols = [c for c in x.columns if not c.startswith("hfd_duration_weeks_")]
assert prev[non_duration_cols].equals(x[non_duration_cols])

# Explicitly guarantee GEO-auto duration evidence is unchanged for all rows.
geo_auto_duration_cols = [
    c for c in x.columns if c.startswith("hfd_duration_weeks_geo_auto_")
]
assert prev[geo_auto_duration_cols].equals(x[geo_auto_duration_cols])

# GSE285614 must remain byte-identical across the entire master row.
assert prev.loc[i614].equals(x.loc[i614])

# GSE324134 must remain byte-identical across the entire master row.
i324 = one("GSE324134")
assert prev.loc[i324].equals(x.loc[i324])

# Previously adjudicated Batch07 target GSE must remain byte-identical.
for gse in [
    "GSE197649", "GSE267123", "GSE274914",
    "GSE283060", "GSE288077", "GSE324134",
]:
    a = prev[prev["GSE"].eq(gse)].sort_values("Contrast_id").reset_index(drop=True)
    b = x[x["GSE"].eq(gse)].sort_values("Contrast_id").reset_index(drop=True)
    assert a.equals(b), f"{gse}: Batch07 canonical row(s) changed unexpectedly"

# ------------------------------------------------------------
# Save Batch08 master
# ------------------------------------------------------------
x.to_csv(OUT, sep="\t", index=False)

# Audit contains both candidates, including the no-change GSE285614 review.
audit = pd.DataFrame([
    {
        "GSE": "GSE285380",
        "Contrast_id": x.loc[i380, "Contrast_id"],
        "prior_final_value": prev.loc[i380, "hfd_duration_weeks_final"],
        "prior_final_confidence": prev.loc[i380, "hfd_duration_weeks_final_confidence"],
        "batch08_final_value": x.loc[i380, "hfd_duration_weeks_final"],
        "batch08_final_confidence": x.loc[i380, "hfd_duration_weeks_final_confidence"],
        "action": "PRIOR_MEDIUM_VALUE_WITHDRAWN_TO_UNRESOLVED",
        "source_citation": citation_380,
        "evidence": evidence_380,
        "note": note_380,
    },
    {
        "GSE": "GSE285614",
        "Contrast_id": x.loc[i614, "Contrast_id"],
        "prior_final_value": prev.loc[i614, "hfd_duration_weeks_final"],
        "prior_final_confidence": prev.loc[i614, "hfd_duration_weeks_final_confidence"],
        "batch08_final_value": x.loc[i614, "hfd_duration_weeks_final"],
        "batch08_final_confidence": x.loc[i614, "hfd_duration_weeks_final_confidence"],
        "action": "MEDIUM_RETAINED_NO_EXACT_WEEK_SOURCE",
        "source_citation": citation_614,
        "evidence": evidence_614,
        "note": (
            "No master-cell change in Batch08. Source directly supports 3 months "
            "for the bulk liver RNA-seq cohort, but no exact week/day value was found; "
            "13.0446-week MEDIUM remains the canonical standardized value."
        ),
    },
])
audit.to_csv(AUDIT_OUT, sep="\t", index=False)

note = f"""HFD duration strict-gap adjudication — Batch08 v1
Canonical input: {MASTER.name}
Output: {OUT.name}

Frozen universe:
- biological eligibility: unchanged
- technical eligibility: unchanged
- GSM membership: unchanged
- contrast membership: unchanged
- sample-role identity: unchanged
- all non-duration moderators: unchanged
- GEO-auto duration evidence: unchanged

GSE285380
- Prior canonical value: 13.0446 wk, MEDIUM.
- GEO GSM treatment metadata report chow/HFHSD for 3 months.
- Publication supplement separately describes a 5-week experimental period followed by sacrifice/liver collection.
- The exact total HFD-start-to-liver-collection exposure for the frozen no-estrogen obese comparator is not explicitly reported.
- Correction: prior 13.0446-wk MEDIUM is withdrawn.
- Final: NOT_REPORTED/UNRESOLVED.
- 3 months + 5 weeks is NOT inferred.

GSE285614
- Publication directly links the bulk liver RNA-seq cohort to 3 months of chow/HFD exposure.
- No exact week/day value is reported.
- Existing standardized 13.0446 wk, MEDIUM is retained.
- No master cell is changed for GSE285614 in Batch08.

GSE324134
- Explicitly untouched in Batch08.
- Canonical Batch07 state remains NOT_REVIEWED.
- GEO-auto duration remains 12 wk, LOW.
- It remains an inclusive-but-not-strict duration candidate, not a reviewed-unresolved case.

Expected coverage after recomputation
- duration inclusive independent GSE: 35
- duration strict independent GSE: 33
- sex strict independent GSE: 37
- strain strict independent GSE: 40
- HFD fat kcal% strict independent GSE: 31
- control diet strict independent GSE: 50

Expected inclusive-but-not-strict duration GSE
- GSE285614: 13.0446 wk, MANUAL_MEDIUM
- GSE324134: 12 wk, GEO_AUTO_LOW

Expected reviewed-unresolved duration GSE
- GSE267072
- GSE285380

Versioning
- Batch07 v1 is never overwritten.
- Batch08 v1 refuses to overwrite an existing Batch08 master.
- Any correction to Batch08 v1 should be emitted as batch08_v1_1 or later.
"""
NOTE_OUT.write_text(note, encoding="utf-8")

print("===== HFD MODERATOR MANUAL ADJUDICATION BATCH08 v1 =====")
print("Canonical input:", MASTER.relative_to(ROOT))
print("Changed contrast:", x.loc[i380, "Contrast_id"])
print("Changed cells:", len(changed))
print("GSE285380:", repr(x.loc[i380, "hfd_duration_weeks_final"]),
      repr(x.loc[i380, "hfd_duration_weeks_final_confidence"]),
      x.loc[i380, "hfd_duration_weeks_adjudication_status"])
print("GSE285614 retained:",
      x.loc[i614, "hfd_duration_weeks_final"],
      x.loc[i614, "hfd_duration_weeks_final_confidence"])
print("GSE324134 untouched:",
      x.loc[i324, "overall_adjudication_status"],
      x.loc[i324, "hfd_duration_weeks_geo_auto_value"],
      x.loc[i324, "hfd_duration_weeks_geo_auto_confidence"])
print("Saved:")
print(" ", OUT.relative_to(ROOT))
print(" ", AUDIT_OUT.relative_to(ROOT))
print(" ", NOTE_OUT.relative_to(ROOT))
print("QA PASS: only GSE285380 duration adjudication/final fields changed.")
