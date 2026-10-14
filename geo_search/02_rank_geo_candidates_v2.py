#!/usr/bin/env python3

import pandas as pd
import re

INFILE = (
    "00_metadata/systematic_search/"
    "geo_candidate_GSE_strict_v2.tsv"
)

AUDITED_FILE = (
    "00_metadata/systematic_search/"
    "already_audited_GSE.txt"
)

OUTFILE = (
    "00_metadata/systematic_search/"
    "strict_screening_ranked_v2.tsv"
)

# ------------------------------------------------------------
# Read input
# ------------------------------------------------------------

df = pd.read_csv(
    INFILE,
    sep="\t",
    header=None,
    names=["GSE", "Title", "Summary"],
    dtype=str
).fillna("")

with open(AUDITED_FILE) as f:
    audited = {
        x.strip()
        for x in f
        if x.strip()
    }


def search(pattern, text):
    return bool(
        re.search(
            pattern,
            text,
            flags=re.I
        )
    )


# ------------------------------------------------------------
# Positive relevance signals
# ------------------------------------------------------------

patterns = {

    "liver":
        r"\b(liver|hepatic)\b",

    "hfd":
        (
            r"\b("
            r"high[- ]fat diet|"
            r"HFD|"
            r"diet[- ]induced obesity|"
            r"DIO"
            r")\b"
        ),

    "control_diet":
        (
            r"\b("
            r"chow|"
            r"normal chow|"
            r"normal diet|"
            r"standard diet|"
            r"control diet|"
            r"low[- ]fat diet|"
            r"LFD|"
            r"NCD|"
            r"regular diet"
            r")\b"
        ),

    "rna_seq":
        (
            r"\b("
            r"RNA[- ]?seq|"
            r"RNA sequencing|"
            r"transcriptom|"
            r"high[- ]throughput sequencing"
            r")\b"
        ),

    "intervention":
        (
            r"\b("
            r"treat|treated|treatment|"
            r"intervention|"
            r"agonist|antagonist|inhibitor|"
            r"drug|compound|"
            r"exercise|training|"
            r"metformin|resveratrol|"
            r"supplement|supplementation"
            r")\b"
        ),

    "recovery":
        (
            r"\b("
            r"weight loss|"
            r"diet switch|"
            r"dietary switch|"
            r"returned? .* chow|"
            r"return .* chow|"
            r"recovery|"
            r"reversal|"
            r"switch(?:ed|ing)? .* diet"
            r")\b"
        ),

    "wt":
        (
            r"\b("
            r"wild[- ]?type|"
            r"\bWT\b|"
            r"C57BL/6|"
            r"C57BL/6J|"
            r"C57BL/6N"
            r")\b"
        )
}


# ------------------------------------------------------------
# Flags requiring later manual separation
# These are NOT automatic exclusions.
# ------------------------------------------------------------

flags = {

    "single_cell":
        (
            r"\b("
            r"single[- ]cell|"
            r"single cell|"
            r"scRNA|"
            r"single[- ]nucleus|"
            r"snRNA"
            r")\b"
        ),

    "isolated_cells":
        (
            r"\b("
            r"isolated hepatocyte|"
            r"primary hepatocyte|"
            r"purified hepatocyte|"
            r"hepatocytes"
            r")\b"
        ),

    "genetic":
        (
            r"\b("
            r"knockout|"
            r"knock-out|"
            r"knock out|"
            r"deficient|"
            r"transgenic|"
            r"conditional knockout|"
            r"KO mice|"
            r"-/-"
            r")\b"
        ),

    "stz":
        r"\b(STZ|streptozotocin)\b",

    "ccl4":
        (
            r"\b("
            r"CCl4|"
            r"carbon tetrachloride"
            r")\b"
        ),

    "special_diet":
        (
            r"\b("
            r"CDAHFD|"
            r"CDAA|"
            r"MCD|"
            r"methionine[- ]choline|"
            r"GAN diet|"
            r"AMLN|"
            r"Western diet|"
            r"high[- ]fat high[- ]sugar|"
            r"HFHS"
            r")\b"
        ),

    "ketogenic":
        r"\b(ketogenic|ketosis|\bKD\b)\b",

    "other_omics":
        (
            r"\b("
            r"ATAC[- ]?seq|"
            r"ChIP[- ]?seq|"
            r"DHS[- ]?seq|"
            r"spatial transcript"
            r")\b"
        )
}


rows = []

for _, r in df.iterrows():

    text = (
        str(r["Title"]) +
        " " +
        str(r["Summary"])
    )

    hit = {
        k: search(p, text)
        for k, p in patterns.items()
    }

    flag = {
        k: search(p, text)
        for k, p in flags.items()
    }

    # --------------------------------------------------------
    # Ranking score
    #
    # This is ONLY manual-review priority.
    # It does not determine inclusion/exclusion.
    # --------------------------------------------------------

    score = 0

    score += 3 if hit["liver"] else 0
    score += 3 if hit["hfd"] else 0
    score += 3 if hit["control_diet"] else 0
    score += 1 if hit["rna_seq"] else 0
    score += 1 if hit["wt"] else 0

    # particularly valuable designs
    score += 2 if hit["intervention"] else 0
    score += 3 if hit["recovery"] else 0

    # lower priority for primary screening,
    # but NEVER automatically exclude
    score -= 2 if flag["single_cell"] else 0
    score -= 1 if flag["isolated_cells"] else 0
    score -= 1 if flag["genetic"] else 0
    score -= 2 if flag["stz"] else 0
    score -= 2 if flag["ccl4"] else 0
    score -= 1 if flag["special_diet"] else 0
    score -= 1 if flag["ketogenic"] else 0

    if score >= 10:
        priority = "HIGH"
    elif score >= 7:
        priority = "MEDIUM"
    else:
        priority = "LOW"

    # --------------------------------------------------------
    # Candidate role
    # --------------------------------------------------------

    if hit["recovery"]:
        role = "RECOVERY_POSSIBLE"

    elif hit["intervention"]:
        role = "INTERVENTION_POSSIBLE"

    elif (
        hit["hfd"] and
        hit["control_diet"]
    ):
        role = "HFD_CORE_POSSIBLE"

    else:
        role = "OTHER_REVIEW"


    rows.append({

        "GSE": r["GSE"],

        "Audit_status":
            (
                "AUDITED"
                if r["GSE"] in audited
                else "NEW"
            ),

        "Priority": priority,
        "Score": score,
        "Candidate_role": role,

        "Liver_hit": hit["liver"],
        "HFD_hit": hit["hfd"],
        "Control_diet_hit":
            hit["control_diet"],
        "WT_hit": hit["wt"],
        "Intervention_hit":
            hit["intervention"],
        "Recovery_hit":
            hit["recovery"],

        "Single_cell_flag":
            flag["single_cell"],
        "Isolated_cell_flag":
            flag["isolated_cells"],
        "Genetic_flag":
            flag["genetic"],
        "STZ_flag":
            flag["stz"],
        "CCl4_flag":
            flag["ccl4"],
        "Special_diet_flag":
            flag["special_diet"],
        "Ketogenic_flag":
            flag["ketogenic"],
        "Other_omics_flag":
            flag["other_omics"],

        "Title": r["Title"],
        "Summary": r["Summary"]
    })


out = pd.DataFrame(rows)

priority_order = {
    "HIGH": 0,
    "MEDIUM": 1,
    "LOW": 2
}

audit_order = {
    "NEW": 0,
    "AUDITED": 1
}

out["_audit"] = (
    out["Audit_status"]
    .map(audit_order)
)

out["_priority"] = (
    out["Priority"]
    .map(priority_order)
)

out = (
    out
    .sort_values(
        [
            "_audit",
            "_priority",
            "Score",
            "GSE"
        ],
        ascending=[
            True,
            True,
            False,
            True
        ]
    )
    .drop(
        columns=[
            "_audit",
            "_priority"
        ]
    )
)


out.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)


print("")
print("========================================")
print("AUTOMATED TRIAGE V2")
print("========================================")

print("")
print("Total GSE:", len(out))

print("")
print("Audit / Priority counts:")
print(
    out.groupby(
        [
            "Audit_status",
            "Priority"
        ]
    ).size()
)

print("")
print("Candidate roles among NEW records:")
print(
    out[
        out["Audit_status"] == "NEW"
    ].groupby(
        "Candidate_role"
    ).size()
    .sort_values(
        ascending=False
    )
)

high = out[
    (out["Audit_status"] == "NEW") &
    (out["Priority"] == "HIGH")
]

print("")
print("NEW HIGH priority:", len(high))
print("")

if len(high) > 0:

    print(
        high[
            [
                "GSE",
                "Score",
                "Candidate_role",
                "Genetic_flag",
                "STZ_flag",
                "Special_diet_flag",
                "Title"
            ]
        ].to_string(
            index=False
        )
    )

print("")
print("Output:")
print(OUTFILE)
