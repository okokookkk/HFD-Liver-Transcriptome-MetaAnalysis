#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd

# ============================================================
# HFD moderator manual adjudication — Batch11 v1
#
# Canonical input:
#   HFD_primary105_moderator_manual_adjudication_batch10_v1.tsv
#
# Scope intentionally restricted to HFD fat kcal%:
#   - GSE229697 -> 60 HIGH
#   - GSE239490 -> 60 HIGH
#
# No other moderator is modified in Batch11, even where the reviewed
# source also exposes duration/strain/control information. This keeps
# the fat-focused pass auditable and prevents cross-field coupling.
#
# Frozen biological/technical eligibility, GSM membership, contrast
# membership, sample-role identity, all GEO-auto fields, and all
# non-fat moderator fields are immutable.
# ============================================================

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=__import__("os").environ.get("HFD_PROJECT_ROOT", "."),
)
args = ap.parse_args()

ROOT = Path(args.root)
OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"

MASTER = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch10_v1.tsv"
OUT = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch11_v1.tsv"
AUDIT_OUT = OUTDIR / "HFD_fat_manual_adjudication_batch11_v1_audit.tsv"
NOTE_OUT = OUTDIR / "HFD_batch11_v1_fat_adjudication_note.txt"

if not MASTER.exists():
    raise FileNotFoundError(MASTER)
if OUT.exists():
    raise FileExistsError(
        f"Refusing to overwrite existing Batch11 master: {OUT}\n"
        "Use batch11_v1_1 or later for any correction."
    )

x = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
prev = x.copy(deep=True)

assert x["GSE"].nunique() == 64
assert x["Contrast_id"].nunique() == 105

VAR = "hfd_fat_kcal_percent"
TARGETS = ["GSE229697", "GSE239490"]

def one(gse):
    z = x[x["GSE"].eq(gse)]
    if len(z) != 1:
        raise RuntimeError(
            f"{gse}: expected exactly one frozen contrast row, observed {len(z)}"
        )
    return z.index[0]

idx = {gse: one(gse) for gse in TARGETS}

def effective(row):
    if row[f"{VAR}_adjudication_status"] != "NOT_REVIEWED":
        return (
            row[f"{VAR}_final"].strip(),
            row[f"{VAR}_final_confidence"].strip(),
            "MANUAL",
        )
    if row[f"{VAR}_geo_auto_status"] == "RESOLVED":
        return (
            row[f"{VAR}_geo_auto_value"].strip(),
            row[f"{VAR}_geo_auto_confidence"].strip(),
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
    val, conf, origin = effective(row)
    if not numeric_ready(val):
        return False
    if origin == "MANUAL":
        return conf == "HIGH"
    if origin == "GEO_AUTO":
        return conf in {"MEDIUM", "HIGH"}
    return False

# Canonical Batch10 anchor.
strict_gses = {
    row["GSE"] for _, row in x.iterrows() if strict_fat(row)
}
assert len(strict_gses) == 39, (
    "Canonical Batch10 fat strict independent-GSE coverage mismatch",
    len(strict_gses), 39
)

# Both Batch11 targets must still be fat pending/missing in canonical Batch10.
for gse, i in idx.items():
    val, conf, origin = effective(x.loc[i])
    assert not numeric_ready(val), (gse, val, conf, origin)
    assert x.loc[i, f"{VAR}_adjudication_status"] == "NOT_REVIEWED", (
        gse, x.loc[i, f"{VAR}_adjudication_status"]
    )

def fill_fat(i, value, source_type, citation, evidence, final_note=""):
    x.loc[i, f"{VAR}_manual_value"] = str(value)
    x.loc[i, f"{VAR}_manual_source_type"] = source_type
    x.loc[i, f"{VAR}_manual_source_citation"] = citation
    x.loc[i, f"{VAR}_manual_evidence"] = evidence
    x.loc[i, f"{VAR}_manual_confidence"] = "HIGH"
    x.loc[i, f"{VAR}_adjudication_status"] = "MANUAL_CONFIRMED"
    x.loc[i, f"{VAR}_final"] = str(value)
    x.loc[i, f"{VAR}_final_confidence"] = "HIGH"
    x.loc[i, f"{VAR}_final_source"] = source_type
    x.loc[i, f"{VAR}_final_note"] = final_note

def set_pub_if_safe(i, pmid, doi, title):
    for col, new in [
        ("publication_pmid", pmid),
        ("publication_doi", doi),
        ("publication_title", title),
    ]:
        old = x.loc[i, col]
        if old and new and old != new:
            raise RuntimeError(
                f"{x.loc[i,'GSE']}: conflicting {col}: "
                f"existing={old!r}, new={new!r}"
            )
        if not old:
            x.loc[i, col] = new

# ============================================================
# GSE229697
#
# Direct HFD GSM evidence:
#   GSM7173369/70/71 treatment protocol:
#   "HFD (60% of energy derived from fat) for 12 weeks"
#   strain C57BL/6N.
#
# Series citation:
#   Lee et al., Clinical and Translational Medicine 2024
#   PMID 39367660; DOI 10.1002/ctm2.1816
#   RNA-seq data explicitly available as GSE229697.
#
# We adjudicate FAT ONLY in Batch11.
# ============================================================
gse = "GSE229697"
i = idx[gse]
citation_229697 = (
    "GEO:GSE229697 selected untreated HFD GSM7173369-GSM7173371; "
    "Lee et al. Clin Transl Med. 2024. PMID:39367660; "
    "DOI:10.1002/ctm2.1816"
)
fill_fat(
    i, "60",
    "GEO_GSM+PUBLICATION",
    citation_229697,
    "Selected untreated HFD liver GSM directly state that C57BL/6N mice "
    "received HFD with 60% of energy derived from fat. Because the source "
    "explicitly reports an energy fraction, the HFD fat moderator is "
    "60 kcal% (energy percentage is invariant to whether energy is expressed "
    "in kcal or kJ). The publication directly links the RNA-seq dataset to "
    "GSE229697.",
    final_note=(
        "Batch11 is fat-only. Direct GSM evidence also contains diet duration, "
        "strain, and normal-chow comparator information, but those moderators "
        "are intentionally not modified in this pass."
    ),
)
set_pub_if_safe(
    i,
    "39367660",
    "10.1002/ctm2.1816",
    "Discovery of a selective cytochrome P450 4A inhibitor for the treatment of metabolic dysfunction-associated fatty liver disease",
)

# ============================================================
# GSE239490
#
# Direct GSM:
#   HFD GSM7666240 etc: "treatment: 16 weeks, 60 %HFD"
#   LFD GSM7666228 etc: "treatment: 16 weeks, 10% LFD"
#
# Primary article Methods:
#   defines HFD = 60 kcal% fat, Research Diets D12492
#   matched LFD = 10 kcal%, Research Diets D12450J
#   mouse RNA-seq accession = GSE239490
#
# The GSM anchors the frozen deposited cohort to the 60% HFD treatment;
# the primary Methods supplies the explicit kcal basis/product definition.
#
# We adjudicate FAT ONLY in Batch11.
# ============================================================
gse = "GSE239490"
i = idx[gse]
citation_239490 = (
    "GEO:GSE239490 selected HFD GSM (e.g. GSM7666240) and LFD GSM "
    "(e.g. GSM7666228); Lutkewitte et al. J Clin Invest. 2024. "
    "PMID:39405118; GEO:GSE239490"
)
fill_fat(
    i, "60",
    "GEO_GSM+PUBLICATION_METHODS",
    citation_239490,
    "Selected deposited HFD liver GSM directly label the treatment as "
    "'16 weeks, 60 %HFD', while matched LFD GSM label '16 weeks, 10% LFD'. "
    "The primary article defines its HFD as Research Diets D12492, "
    "60 kcal% fat, and explicitly states that the mouse RNA-seq data were "
    "deposited as GSE239490. Together these provide a direct cohort anchor "
    "and explicit caloric-fat definition for final HFD fat = 60 kcal%.",
    final_note=(
        "Batch11 is fat-only. The 16-week duration and matched LFD are present "
        "in deposited GSM metadata but are intentionally not modified here."
    ),
)
# DOI intentionally not hard-coded because PMID/GSE linkage is sufficient
# and avoids introducing an unchecked bibliographic field.
set_pub_if_safe(
    i,
    "39405118",
    "",
    "Adipocyte lipin 1 expression associates with human metabolic health and regulates systemic metabolism in mice",
)

# ------------------------------------------------------------
# Immutability QA
# ------------------------------------------------------------

# All GEO-auto fields must remain byte-identical.
geo_auto_cols = [c for c in x.columns if "_geo_auto_" in c]
assert prev[geo_auto_cols].equals(x[geo_auto_cols])

# Frozen identity/membership fields immutable.
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

# All non-fat moderator columns must remain identical for all rows.
other_mod_prefixes = [
    "hfd_duration_weeks_",
    "sex_",
    "strain_",
    "control_diet_type_",
    "collection_ZT_",
]
for prefix in other_mod_prefixes:
    cols = [c for c in x.columns if c.startswith(prefix)]
    if cols:
        assert prev[cols].equals(x[cols]), prefix

# Only target fat adjudication/final columns and publication metadata may change.
allowed = {
    c for c in x.columns
    if c.startswith(VAR + "_") and "_geo_auto_" not in c
} | {"publication_pmid", "publication_doi", "publication_title"}

changed = []
for row_i in x.index:
    gse = x.loc[row_i, "GSE"]
    for col in x.columns:
        if prev.loc[row_i, col] != x.loc[row_i, col]:
            changed.append((gse, x.loc[row_i, "Contrast_id"], col))
            assert gse in TARGETS, ("Unexpected GSE modified", gse, col)
            assert col in allowed, ("Unexpected field modified", gse, col)

assert {g for g, _, _ in changed} == set(TARGETS)

# Untargeted GSE rows must remain fully byte-identical.
untargeted = ~x["GSE"].isin(TARGETS)
assert prev.loc[untargeted].reset_index(drop=True).equals(
    x.loc[untargeted].reset_index(drop=True)
)

# Final target QA.
for gse in TARGETS:
    i = idx[gse]
    assert x.loc[i, f"{VAR}_final"] == "60"
    assert x.loc[i, f"{VAR}_final_confidence"] == "HIGH"
    assert x.loc[i, f"{VAR}_adjudication_status"] == "MANUAL_CONFIRMED"

x.to_csv(OUT, sep="\t", index=False)

audit = pd.DataFrame([
    {
        "GSE": gse,
        "Contrast_id": x.loc[idx[gse], "Contrast_id"],
        "Moderator": VAR,
        "Prior_final_value": prev.loc[idx[gse], f"{VAR}_final"],
        "Prior_final_confidence": prev.loc[idx[gse], f"{VAR}_final_confidence"],
        "Batch11_final_value": x.loc[idx[gse], f"{VAR}_final"],
        "Batch11_final_confidence": x.loc[idx[gse], f"{VAR}_final_confidence"],
        "Source_type": x.loc[idx[gse], f"{VAR}_manual_source_type"],
        "Source_citation": x.loc[idx[gse], f"{VAR}_manual_source_citation"],
        "Evidence": x.loc[idx[gse], f"{VAR}_manual_evidence"],
        "Final_note": x.loc[idx[gse], f"{VAR}_final_note"],
        "Action": "MANUAL_HIGH_DIRECT_SOURCE",
    }
    for gse in TARGETS
])
audit.to_csv(AUDIT_OUT, sep="\t", index=False)

NOTE_OUT.write_text(
f"""HFD moderator manual adjudication — Batch11 v1
Canonical input: {MASTER.name}
Output: {OUT.name}

Scope
- HFD fat kcal% only.
- No duration, sex, strain, control-diet, sample-role, eligibility,
  GSM membership, or contrast membership changes.

Direct HIGH fat upgrades
- GSE229697: 60 kcal% HIGH.
  Untreated HFD GSM directly state 60% of energy derived from fat.
  Publication directly links RNA-seq to GSE229697.
- GSE239490: 60 kcal% HIGH.
  Deposited HFD GSM directly label the cohort 60% HFD; the primary
  publication defines HFD as D12492, 60 kcal% fat and directly links
  mouse RNA-seq to GSE239490.

Deliberately deferred despite direct ancillary metadata
- GSE229697: duration/strain/control diet not changed in this fat-only pass.
- GSE239490: duration/control diet not changed in this fat-only pass.

Expected independent-GSE coverage after recomputation
- HFD duration strict: 40 (unchanged)
- sex strict: 42 (unchanged)
- strain strict: 41 (unchanged)
- HFD fat kcal% inclusive: 42
- HFD fat kcal% strict: 41
- control diet strict: 52 (unchanged)

Versioning
- Batch10 v1 is not overwritten.
- Batch11 v1 refuses overwrite.
- Corrections must supersede as Batch11 v1.1 or later.
""",
    encoding="utf-8",
)

print("===== HFD MODERATOR MANUAL ADJUDICATION BATCH11 v1 =====")
print("Canonical input:", MASTER.relative_to(ROOT))
print("Target GSE:", ", ".join(TARGETS))
print("Changed cells:", len(changed))
print()
for gse in TARGETS:
    i = idx[gse]
    print(
        gse,
        "fat=",
        x.loc[i, f"{VAR}_final"],
        x.loc[i, f"{VAR}_final_confidence"],
    )
print("\nSaved:")
print(" ", OUT.relative_to(ROOT))
print(" ", AUDIT_OUT.relative_to(ROOT))
print(" ", NOTE_OUT.relative_to(ROOT))
print("QA PASS: fat-only target changes; all other moderator/frozen fields unchanged.")
