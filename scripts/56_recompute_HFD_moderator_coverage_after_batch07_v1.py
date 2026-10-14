#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

MODS = ["hfd_duration_weeks", "sex", "strain", "hfd_fat_kcal_percent", "control_diet_type"]
VALID_STRAINS = {
    "C57BL/6J", "C57BL/6", "C57BL/6NHsd", "129S1/SvImJ",
    "C57BL/6N", "FVB.129P2", "C57BL/6NTac"
}
VALID_SEX = {"MALE", "FEMALE"}
VALID_CONTROL = {"CHOW", "LFD", "NORMAL_DIET", "CONTROL_DIET_UNSPECIFIED"}

def effective(r, mod):
    if r[f"{mod}_adjudication_status"] != "NOT_REVIEWED":
        return r[f"{mod}_final"].strip(), r[f"{mod}_final_confidence"].strip(), "MANUAL"
    if r[f"{mod}_geo_auto_status"] == "RESOLVED":
        return r[f"{mod}_geo_auto_value"].strip(), r[f"{mod}_geo_auto_confidence"].strip(), "GEO_AUTO"
    return "", "", "NONE"

def ready(mod, val):
    if not val or val in {"UNKNOWN", "UNRESOLVED", "NOT_REPORTED", "NA", "N/A"}:
        return False
    if mod == "hfd_duration_weeks":
        try:
            return float(val) > 0
        except ValueError:
            return False
    if mod == "hfd_fat_kcal_percent":
        try:
            x = float(val)
            return 0 < x <= 100
        except ValueError:
            return False
    if mod == "sex":
        return val in VALID_SEX
    if mod == "strain":
        return val in VALID_STRAINS
    if mod == "control_diet_type":
        return val in VALID_CONTROL
    return False

def strict_ready(mod, val, conf, origin):
    if not ready(mod, val):
        return False
    if origin == "MANUAL":
        return conf == "HIGH"
    if origin == "GEO_AUTO":
        return conf in {"MEDIUM", "HIGH"}
    return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--master",
        default="00_metadata/contrasts/moderators/manual_adjudication/"
                "HFD_primary105_moderator_manual_adjudication_batch07_v1.tsv",
    )
    ap.add_argument(
        "--out",
        default="00_metadata/contrasts/moderators/manual_adjudication/"
                "HFD_primary105_moderator_coverage_after_batch07_v1.tsv",
    )
    args = ap.parse_args()

    with open(args.master, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    out_rows = []
    for mod in MODS:
        known = inclusive = strict = unresolved = pending = 0
        inclusive_gse, strict_gse = set(), set()
        for r in rows:
            val, conf, origin = effective(r, mod)
            is_ready = ready(mod, val)
            is_strict = strict_ready(mod, val, conf, origin)
            if val:
                known += 1
            if is_ready:
                inclusive += 1
                inclusive_gse.add(r["GSE"])
            if is_strict:
                strict += 1
                strict_gse.add(r["GSE"])
            if r[f"{mod}_manual_confidence"] == "UNRESOLVED":
                unresolved += 1
            if r[f"{mod}_adjudication_status"] == "NOT_REVIEWED" and not is_ready:
                pending += 1
        out_rows.append({
            "Moderator": mod,
            "Total_contrasts": len(rows),
            "Known_contrasts": known,
            "Inclusive_ready_contrasts": inclusive,
            "Independent_GSE_inclusive_ready": len(inclusive_gse),
            "Strict_ready_contrasts": strict,
            "Independent_GSE_strict_ready": len(strict_gse),
            "Reviewed_but_unresolved_contrasts": unresolved,
            "Pending_review_contrasts": pending,
        })

    fields = [
        "Moderator","Total_contrasts","Known_contrasts","Inclusive_ready_contrasts",
        "Independent_GSE_inclusive_ready","Strict_ready_contrasts",
        "Independent_GSE_strict_ready","Reviewed_but_unresolved_contrasts",
        "Pending_review_contrasts"
    ]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)

    print("===== HFD MODERATOR COVERAGE AFTER BATCH07 v1 =====")
    for x in out_rows:
        print(
            f'{x["Moderator"]:25s} total={x["Total_contrasts"]:3d} '
            f'known={x["Known_contrasts"]:3d} '
            f'inclusive={x["Inclusive_ready_contrasts"]:3d}/GSE{x["Independent_GSE_inclusive_ready"]:2d} '
            f'strict={x["Strict_ready_contrasts"]:3d}/GSE{x["Independent_GSE_strict_ready"]:2d} '
            f'unresolved={x["Reviewed_but_unresolved_contrasts"]:2d} '
            f'pending={x["Pending_review_contrasts"]:2d}'
        )
    print(f"Wrote: {out}")

if __name__ == "__main__":
    main()
