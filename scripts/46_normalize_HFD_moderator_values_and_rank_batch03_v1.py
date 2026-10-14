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
INFILE = (
    ROOT /
    "00_metadata/contrasts/moderators/manual_adjudication/"
    "HFD_primary105_moderator_effective_after_batch02_v1.tsv"
)
MASTER = (
    ROOT /
    "00_metadata/contrasts/moderators/manual_adjudication/"
    "HFD_primary105_moderator_manual_adjudication_batch02_v1.tsv"
)
OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"

if not INFILE.exists():
    raise FileNotFoundError(INFILE)
if not MASTER.exists():
    raise FileNotFoundError(MASTER)

e = pd.read_csv(INFILE, sep="\t", dtype=str).fillna("")
m = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")

assert e["Contrast_id"].nunique() == 105
assert e["GSE"].nunique() == 64

vars_ = [
    "hfd_duration_weeks",
    "sex",
    "strain",
    "hfd_fat_kcal_percent",
    "control_diet_type",
]

def canon_number(v):
    v = str(v).strip()
    if v in {"", "UNKNOWN", "CONFLICT", "MIXED"}:
        return v
    try:
        f = float(v)
    except Exception:
        return v
    if f.is_integer():
        return str(int(f))
    return f"{f:g}"

# ------------------------------------------------------------
# Canonicalize numeric moderators only.
# Preserve all source/raw/final fields in master unchanged.
# ------------------------------------------------------------
e2 = e.copy()
for var in ["hfd_duration_weeks", "hfd_fat_kcal_percent"]:
    col = f"{var}_effective"
    e2[col] = e2[col].map(canon_number)

# ------------------------------------------------------------
# Corrected coverage/distribution
# ------------------------------------------------------------
coverage_rows = []
distribution_rows = []

for var in vars_:
    known = int(e2[f"{var}_known"].astype(int).sum())
    ready = int(e2[f"{var}_analysis_ready"].astype(int).sum())
    ready_gse = e2.loc[
        e2[f"{var}_analysis_ready"].eq("1"), "GSE"
    ].nunique()

    coverage_rows.append({
        "Moderator": var,
        "Total_contrasts": len(e2),
        "Known_contrasts": known,
        "Known_percent": 100 * known / len(e2),
        "Analysis_ready_contrasts": ready,
        "Analysis_ready_percent": 100 * ready / len(e2),
        "Independent_GSE_with_ready_value": ready_gse,
    })

    z = e2[e2[f"{var}_known"].eq("1")].copy()
    for value, q in z.groupby(f"{var}_effective", sort=True):
        distribution_rows.append({
            "Moderator": var,
            "Value": value,
            "Contrasts": q["Contrast_id"].nunique(),
            "Independent_GSE": q["GSE"].nunique(),
            "Analysis_ready_contrasts": int(
                q[f"{var}_analysis_ready"].astype(int).sum()
            ),
        })

coverage = pd.DataFrame(coverage_rows)
dist = pd.DataFrame(distribution_rows)

# ------------------------------------------------------------
# Independent-GSE marginal-gain curation queue.
#
# Goal: do not overvalue a 12-contrast GSE merely because it has many
# dependent contrasts. Score each unresolved moderator once per GSE.
#
# Priority weights reflect current scientific plan:
# duration > sex > strain > HFD fat > control diet.
# ------------------------------------------------------------
weights = {
    "hfd_duration_weeks": 5,
    "sex": 4,
    "strain": 3,
    "hfd_fat_kcal_percent": 2,
    "control_diet_type": 1,
}

queue_rows = []

for gse, z in e2.groupby("GSE", sort=True):
    row = {
        "GSE": gse,
        "n_contrasts": z["Contrast_id"].nunique(),
    }

    score = 0
    unresolved_vars = []
    ready_vars = []

    for var in vars_:
        # GSE gets credit as "ready" if at least one contrast has a ready value.
        # It needs review if any contrast is not ready, because multi-stratum
        # GSEs may require contrast-specific adjudication.
        n_not_ready = int((z[f"{var}_analysis_ready"] != "1").sum())
        n_ready = int((z[f"{var}_analysis_ready"] == "1").sum())

        row[f"{var}_not_ready_n"] = n_not_ready
        row[f"{var}_ready_n"] = n_ready

        if n_not_ready > 0:
            score += weights[var]
            unresolved_vars.append(var)
        else:
            ready_vars.append(var)

    row["gse_marginal_gain_score"] = score
    row["unresolved_moderators"] = ",".join(unresolved_vars)
    row["fully_ready_moderators"] = ",".join(ready_vars)

    # Prefer papers that can resolve several moderator dimensions at once.
    row["n_unresolved_moderator_types"] = len(unresolved_vars)

    # Existing manual status.
    mz = m[m["GSE"].eq(gse)]
    row["already_manual_adjudicated"] = int(
        len(mz) > 0 and
        mz["overall_adjudication_status"]
        .str.contains("ADJUDICATED", na=False)
        .all()
    )

    queue_rows.append(row)

queue = pd.DataFrame(queue_rows)
queue = queue[
    (queue["gse_marginal_gain_score"] > 0)
    & (queue["already_manual_adjudicated"] == 0)
].copy()

queue = queue.sort_values(
    [
        "gse_marginal_gain_score",
        "n_unresolved_moderator_types",
        "n_contrasts",
        "GSE",
    ],
    ascending=[False, False, True, True],
)

# Recommended batch03:
# Favor high marginal gain and independent GSE diversity.
batch03 = queue.head(12).copy()
batch03.insert(0, "Batch03_order", range(1, len(batch03) + 1))

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
eff_out = OUTDIR / "HFD_primary105_moderator_effective_after_batch02_v1_1.tsv"
cov_out = OUTDIR / "HFD_primary105_moderator_coverage_after_batch02_v1_1.tsv"
dist_out = OUTDIR / "HFD_primary105_moderator_distribution_after_batch02_v1_1.tsv"
queue_out = OUTDIR / "HFD_primary64_GSE_marginal_gain_curation_queue_after_batch02_v1.tsv"
batch_out = OUTDIR / "HFD_batch03_recommended_GSE_queue_v1.tsv"

e2.to_csv(eff_out, sep="\t", index=False)
coverage.to_csv(cov_out, sep="\t", index=False)
dist.to_csv(dist_out, sep="\t", index=False)
queue.to_csv(queue_out, sep="\t", index=False)
batch03.to_csv(batch_out, sep="\t", index=False)

print("===== HFD MODERATOR COVERAGE AFTER BATCH02 v1.1 =====")
print("\n===== EFFECTIVE COVERAGE =====")
print(coverage.to_string(index=False))

print("\n===== CANONICALIZED KEY DISTRIBUTIONS =====")
for var in vars_:
    print(f"\n--- {var} ---")
    z = dist[dist["Moderator"].eq(var)]
    print(
        z.sort_values(
            ["Independent_GSE", "Contrasts"],
            ascending=False
        ).to_string(index=False)
    )

print("\n===== TOP GSE BY INDEPENDENT MARGINAL GAIN =====")
cols = [
    "GSE",
    "n_contrasts",
    "gse_marginal_gain_score",
    "n_unresolved_moderator_types",
    "unresolved_moderators",
]
print(queue[cols].head(20).to_string(index=False))

print("\n===== RECOMMENDED BATCH03 =====")
print(batch03[["Batch03_order"] + cols].to_string(index=False))

print("\nSaved:")
for p in [eff_out, cov_out, dist_out, queue_out, batch_out]:
    print(p.relative_to(ROOT))

print("\nNOTE:")
print("- Numeric canonicalization changes only reporting categories (e.g. 12.0 -> 12).")
print("- No manual adjudication, biological label, contrast membership, or eligibility changed.")
print("- Batch03 queue is scored once per unresolved moderator per GSE, not per contrast.")
