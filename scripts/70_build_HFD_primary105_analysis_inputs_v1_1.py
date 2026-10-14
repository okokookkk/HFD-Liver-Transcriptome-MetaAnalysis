#!/usr/bin/env python3
from pathlib import Path
import hashlib
import pandas as pd

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

MANIFEST = ROOT / "00_metadata/technical_harmonization/frozen_v1/HFD_primary_ENA_R1R2_download_manifest_frozen_v1.tsv"
FINAL_STATUS = ROOT / "00_metadata/processing_production/final_audit_v1/HFD_primary64_final_resolved_status_v1.tsv"
MODERATOR_MASTER = ROOT / "00_metadata/contrasts/moderators/manual_adjudication/HFD_primary105_moderator_manual_adjudication_batch11_v1.tsv"

OUTDIR = ROOT / "00_metadata/analysis_inputs/frozen_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_MEMBERSHIP = OUTDIR / "HFD_primary105_contrast_GSM_membership_v1_1.tsv"
OUT_CONTRAST = OUTDIR / "HFD_primary105_contrast_analysis_input_audit_v1_1.tsv"
OUT_GSE = OUTDIR / "HFD_primary64_GSE_contrast_multiplicity_v1_1.tsv"
OUT_SHARED = OUTDIR / "HFD_primary105_shared_GSM_audit_v1_1.tsv"
OUT_MOD = OUTDIR / "HFD_primary105_moderator_join_v1_1.tsv"
OUT_SUMMARY = OUTDIR / "HFD_primary105_analysis_input_audit_summary_v1_1.txt"

for q in [MANIFEST, FINAL_STATUS, MODERATOR_MASTER]:
    if not q.exists():
        raise RuntimeError(f"Required input not found: {q}")

m = pd.read_csv(MANIFEST, sep="\t", dtype=str).fillna("")
fs = pd.read_csv(FINAL_STATUS, sep="\t", dtype=str).fillna("")
mod = pd.read_csv(MODERATOR_MASTER, sep="\t", dtype=str).fillna("")

required_manifest_cols = {
    "GSE", "GSM", "SRR", "Contrast_id", "Role",
    "FASTQ_component", "FASTQ_filename"
}
missing = required_manifest_cols - set(m.columns)
if missing:
    raise RuntimeError(f"Frozen manifest missing columns: {sorted(missing)}")

expected = {
    "GSE": 64,
    "GSM": 837,
    "SRR": 951,
    "Contrast_id": 105,
}
for col, n in expected.items():
    observed = m[col].nunique()
    if observed != n:
        raise RuntimeError(f"Expected {n} unique {col}, observed {observed}")

if len(m) != 1902:
    raise RuntimeError(f"Expected 1902 FASTQ rows, observed {len(m)}")

if set(m["Role"].unique()) != {"HFD", "CONTROL"}:
    raise RuntimeError(f"Unexpected Role values: {sorted(m['Role'].unique())}")

if fs["GSE"].nunique() != 64 or len(fs) != 64:
    raise RuntimeError("Final resolved production status must contain exactly 64 unique GSE")
if not (fs["Final_status"] == "PASS").all():
    raise RuntimeError("Final production audit contains non-PASS GSE")

membership = (
    m[["GSE", "Contrast_id", "Role", "GSM", "SRR"]]
    .drop_duplicates()
    .sort_values(["GSE", "Contrast_id", "Role", "GSM", "SRR"])
    .reset_index(drop=True)
)

contrast_gse_n = membership.groupby("Contrast_id")["GSE"].nunique()
if not (contrast_gse_n == 1).all():
    raise RuntimeError(
        "At least one Contrast_id maps to multiple GSE:\n"
        + contrast_gse_n[contrast_gse_n != 1].to_string()
    )

# Reuse is audited rather than assumed absent.
gsm_contrast_n = (
    membership[["GSM", "Contrast_id"]]
    .drop_duplicates()
    .groupby("GSM")["Contrast_id"]
    .nunique()
    .rename("N_contrasts_using_GSM")
)

gsm_role_n = (
    membership[["GSM", "Role"]]
    .drop_duplicates()
    .groupby("GSM")["Role"]
    .nunique()
    .rename("N_roles_for_GSM")
)

role_conflict = gsm_role_n[gsm_role_n > 1]
if len(role_conflict):
    raise RuntimeError(
        "At least one GSM changes HFD/CONTROL role across contrasts:\n"
        + role_conflict.to_string()
    )

membership2 = (
    membership[["GSE", "Contrast_id", "Role", "GSM"]]
    .drop_duplicates()
    .merge(gsm_contrast_n, on="GSM", how="left")
    .sort_values(["GSE", "Contrast_id", "Role", "GSM"])
)

membership2.to_csv(OUT_MEMBERSHIP, sep="\t", index=False)

shared = membership2[membership2["N_contrasts_using_GSM"] > 1].copy()
shared.to_csv(OUT_SHARED, sep="\t", index=False)

contrast_rows = []

for contrast_id, z in membership2.groupby("Contrast_id", sort=True):
    gse_values = z["GSE"].unique()
    if len(gse_values) != 1:
        raise RuntimeError(f"{contrast_id}: non-unique GSE")
    gse = gse_values[0]

    hfd_gsm = sorted(z.loc[z["Role"] == "HFD", "GSM"].unique())
    ctl_gsm = sorted(z.loc[z["Role"] == "CONTROL", "GSM"].unique())

    if len(hfd_gsm) < 1 or len(ctl_gsm) < 1:
        raise RuntimeError(
            f"{contrast_id}: HFD={len(hfd_gsm)} CONTROL={len(ctl_gsm)}"
        )

    matrix = ROOT / "03_counts/gsm" / gse / f"{gse}_GSM_counts_v1.tsv"
    if not matrix.exists():
        raise RuntimeError(f"{contrast_id}: missing final matrix {matrix}")

    cols = pd.read_csv(matrix, sep="\t", nrows=0).columns.tolist()
    if "Geneid" not in cols:
        raise RuntimeError(f"{gse}: matrix lacks Geneid")
    observed_gsm = set(cols) - {"Geneid"}

    required_gsm = set(hfd_gsm) | set(ctl_gsm)
    missing_gsm = sorted(required_gsm - observed_gsm)
    if missing_gsm:
        raise RuntimeError(
            f"{contrast_id}: frozen GSM missing from final matrix: {missing_gsm}"
        )

    sr = membership[membership["Contrast_id"].eq(contrast_id)]
    hfd_srr = sr.loc[sr["Role"] == "HFD", "SRR"].nunique()
    ctl_srr = sr.loc[sr["Role"] == "CONTROL", "SRR"].nunique()

    hfd_shared = int(
        (gsm_contrast_n.reindex(hfd_gsm).fillna(0).astype(int) > 1).sum()
    )
    ctl_shared = int(
        (gsm_contrast_n.reindex(ctl_gsm).fillna(0).astype(int) > 1).sum()
    )

    contrast_rows.append({
        "Contrast_id": contrast_id,
        "GSE": gse,
        "N_HFD_GSM": len(hfd_gsm),
        "N_CONTROL_GSM": len(ctl_gsm),
        "N_TOTAL_GSM": len(required_gsm),
        "N_HFD_SRR": hfd_srr,
        "N_CONTROL_SRR": ctl_srr,
        "N_TOTAL_SRR": hfd_srr + ctl_srr,
        "N_HFD_GSM_shared_across_contrasts": hfd_shared,
        "N_CONTROL_GSM_shared_across_contrasts": ctl_shared,
        "Any_shared_GSM": "YES" if (hfd_shared + ctl_shared) > 0 else "NO",
        "HFD_GSMs": ",".join(hfd_gsm),
        "CONTROL_GSMs": ",".join(ctl_gsm),
        "Final_matrix": str(matrix.relative_to(ROOT)),
        "Constructable": "PASS",
    })

contrast = pd.DataFrame(contrast_rows).sort_values(["GSE", "Contrast_id"])

if len(contrast) != 105 or contrast["Contrast_id"].nunique() != 105:
    raise RuntimeError("Did not reconstruct exactly 105 contrasts")

gse_rows = []
for gse, z in contrast.groupby("GSE", sort=True):
    cm = membership2[membership2["GSE"].eq(gse)]
    ids = sorted(z["Contrast_id"].tolist())
    gse_rows.append({
        "GSE": gse,
        "N_contrasts": len(ids),
        "N_unique_GSM_used": cm["GSM"].nunique(),
        "N_HFD_GSM_memberships": (
            cm.loc[cm["Role"] == "HFD", ["Contrast_id", "GSM"]]
            .drop_duplicates()
            .shape[0]
        ),
        "N_CONTROL_GSM_memberships": (
            cm.loc[cm["Role"] == "CONTROL", ["Contrast_id", "GSM"]]
            .drop_duplicates()
            .shape[0]
        ),
        "N_shared_GSM_within_GSE": cm.loc[
            cm["N_contrasts_using_GSM"] > 1, "GSM"
        ].nunique(),
        "Contrast_ids": ",".join(ids),
    })

gse_mult = pd.DataFrame(gse_rows).sort_values(
    ["N_contrasts", "GSE"], ascending=[False, True]
)

if len(gse_mult) != 64:
    raise RuntimeError(f"Expected 64 GSE multiplicity rows, observed {len(gse_mult)}")

gse_mult.to_csv(OUT_GSE, sep="\t", index=False)

if "Contrast_id" not in mod.columns:
    raise RuntimeError(
        f"Moderator master lacks Contrast_id. Columns={list(mod.columns)}"
    )

if len(mod) != 105 or mod["Contrast_id"].nunique() != 105:
    raise RuntimeError(
        f"Moderator master expected 105 unique rows; "
        f"rows={len(mod)} unique={mod['Contrast_id'].nunique()}"
    )

expected_mod_cols = [
    "hfd_duration_weeks_final",
    "sex_final",
    "strain_final",
    "hfd_fat_kcal_percent_final",
    "control_diet_type_final",
]
missing_mod_cols = [c for c in expected_mod_cols if c not in mod.columns]
if missing_mod_cols:
    raise RuntimeError(
        f"Moderator master missing expected columns: {missing_mod_cols}"
    )

contrast_ids = set(contrast["Contrast_id"])
mod_ids = set(mod["Contrast_id"])
if contrast_ids != mod_ids:
    raise RuntimeError(
        "Moderator Contrast_id set differs from frozen contrasts.\n"
        f"Missing={sorted(contrast_ids - mod_ids)}\n"
        f"Extra={sorted(mod_ids - contrast_ids)}"
    )

joined = contrast.merge(
    mod,
    on="Contrast_id",
    how="left",
    validate="one_to_one",
    suffixes=("", "_moderator"),
)

# Canonical analysis aliases. These are copies of the adjudicated *_final
# columns; the original 108-column moderator audit record is retained.
joined["hfd_duration_weeks"] = joined["hfd_duration_weeks_final"]
joined["sex"] = joined["sex_final"]
joined["strain"] = joined["strain_final"]
joined["hfd_fat_kcal_percent"] = joined["hfd_fat_kcal_percent_final"]
joined["control_diet_type"] = joined["control_diet_type_final"]

# Preserve collection time as an optional sixth moderator/covariate when available.
if "collection_ZT_final" in joined.columns:
    joined["collection_ZT"] = joined["collection_ZT_final"]

if "GSE_moderator" in joined.columns:
    mismatch = joined[
        joined["GSE_moderator"].ne("")
        & joined["GSE"].ne(joined["GSE_moderator"])
    ]
    if len(mismatch):
        raise RuntimeError(
            "Moderator GSE disagrees with frozen GSE:\n"
            + mismatch[
                ["Contrast_id", "GSE", "GSE_moderator"]
            ].to_string(index=False)
        )

joined.to_csv(OUT_MOD, sep="\t", index=False)
joined.to_csv(OUT_CONTRAST, sep="\t", index=False)

shared_gsm_n = int((gsm_contrast_n > 1).sum())
max_reuse = int(gsm_contrast_n.max())
gse_multi_n = int((gse_mult["N_contrasts"] > 1).sum())
max_contrasts_per_gse = int(gse_mult["N_contrasts"].max())

def known_n(col):
    return int((joined[col].astype(str).str.strip() != "").sum())

membership_hash = hashlib.sha256(OUT_MEMBERSHIP.read_bytes()).hexdigest()

summary = [
    "===== HFD PRIMARY105 ANALYSIS-INPUT AUDIT v1.1 =====",
    "",
    "FINAL_STATUS=PASS",
    "Production_GSE=64",
    "Frozen_contrasts=105",
    "Frozen_GSM=837",
    "Frozen_SRR=951",
    "Frozen_FASTQ=1902",
    "Constructable_contrasts=105/105",
    "",
    f"GSE_with_multiple_contrasts={gse_multi_n}/64",
    f"Maximum_contrasts_per_GSE={max_contrasts_per_gse}",
    f"GSM_reused_across_contrasts={shared_gsm_n}",
    f"Maximum_contrasts_using_one_GSM={max_reuse}",
    f"Contrasts_with_any_shared_GSM={int((contrast['Any_shared_GSM'] == 'YES').sum())}/105",
    "",
    "Moderator_master=HFD_primary105_moderator_manual_adjudication_batch11_v1.tsv",
    "Moderator_Contrast_id_join=105/105 PASS",
    f"Nonblank_hfd_duration_weeks_final={known_n('hfd_duration_weeks_final')}/105",
    f"Nonblank_sex_final={known_n('sex_final')}/105",
    f"Nonblank_strain_final={known_n('strain_final')}/105",
    f"Nonblank_hfd_fat_kcal_percent_final={known_n('hfd_fat_kcal_percent_final')}/105",
    f"Nonblank_control_diet_type_final={known_n('control_diet_type_final')}/105",
    "",
    f"Membership_SHA256={membership_hash}",
    "",
    f"Contrast_membership={OUT_MEMBERSHIP.relative_to(ROOT)}",
    f"Contrast_analysis_audit={OUT_CONTRAST.relative_to(ROOT)}",
    f"GSE_contrast_multiplicity={OUT_GSE.relative_to(ROOT)}",
    f"Shared_GSM_audit={OUT_SHARED.relative_to(ROOT)}",
    f"Moderator_join={OUT_MOD.relative_to(ROOT)}",
    "",
    "NO_DIFFERENTIAL_EXPRESSION_PERFORMED=YES",
    "FROZEN_MEMBERSHIP_MODIFIED=NO",
    "MODERATOR_MASTER_MODIFIED=NO",
    "FINAL_COUNT_MATRICES_MODIFIED=NO",
]

OUT_SUMMARY.write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))

