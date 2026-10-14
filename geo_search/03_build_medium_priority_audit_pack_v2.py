#!/usr/bin/env python3

import pandas as pd
import subprocess
import pathlib
import time
import re
import sys

BASE = pathlib.Path(
    "00_metadata/systematic_search"
)

RANKED = BASE / "strict_screening_ranked_v2.tsv"
SERIES_META = BASE / "geo_candidate_GSE_strict_v2.tsv"

CACHE = BASE / "audit_cache"
CACHE.mkdir(
    parents=True,
    exist_ok=True
)

OUT_SERIES = BASE / "medium_priority_series_v2.tsv"
OUT_SAMPLES = BASE / "medium_priority_samples_v2.tsv"
OUT_DESIGN = BASE / "medium_priority_design_summary_v2.tsv"
OUT_FAILED = BASE / "medium_priority_failed_GSE.txt"


# ============================================================
# 1. Read ranked candidates
# ============================================================

ranked = pd.read_csv(
    RANKED,
    sep="\t",
    dtype=str
).fillna("")


targets = (
    ranked[
        (ranked["Audit_status"] == "NEW") &
        (ranked["Priority"] == "MEDIUM")
    ]["GSE"]
    .drop_duplicates()
    .tolist()
)


for gse in targets:

    if not re.fullmatch(
        r"GSE[0-9]+",
        gse
    ):
        raise ValueError(
            f"Invalid accession: {gse}"
        )


print(
    "========================================"
)
print(
    "MEDIUM-PRIORITY AUDIT PACK V2"
)
print(
    "========================================"
)
print(
    f"Target GSE: {len(targets)}"
)
print()


# ============================================================
# 2. Read local Series metadata
# ============================================================

series_meta = pd.read_csv(
    SERIES_META,
    sep="\t",
    header=None,
    names=[
        "GSE",
        "Series_title",
        "Series_summary"
    ],
    dtype=str
).fillna("")


series_meta = (
    series_meta
    .drop_duplicates(
        subset="GSE"
    )
    .set_index("GSE")
)


# ============================================================
# 3. Download / cache GEO summaries
# ============================================================

failed = []


def cache_is_valid(
    cache_file,
    gse
):

    if not cache_file.exists():
        return False

    try:

        with open(
            cache_file,
            encoding="utf-8",
            errors="replace"
        ) as f:

            for line in f:

                first = (
                    line.rstrip("\n")
                    .split("\t")[0]
                )

                if first == gse:
                    return True

    except Exception:
        return False

    return False


for i, gse in enumerate(
    targets,
    start=1
):

    cache_file = (
        CACHE /
        f"{gse}.tsv"
    )


    if cache_is_valid(
        cache_file,
        gse
    ):

        print(
            f"[{i:02d}/{len(targets)}] "
            f"{gse} : cached"
        )

        continue


    print(
        f"[{i:02d}/{len(targets)}] "
        f"{gse} : retrieving"
    )


    command = f"""
set -o pipefail

esearch \
    -db gds \
    -query '{gse}[ACCN]' \
| esummary \
| xtract \
    -pattern DocumentSummary \
    -element \
        DocumentSummary/Accession \
        DocumentSummary/title \
        DocumentSummary/summary
"""


    success = False


    for attempt in range(
        1,
        7
    ):

        try:

            result = subprocess.run(
                [
                    "bash",
                    "-c",
                    command
                ],
                capture_output=True,
                text=True,
                timeout=180
            )

        except subprocess.TimeoutExpired:

            print(
                f"    attempt {attempt}: timeout"
            )

            time.sleep(
                attempt * 3
            )

            continue


        if result.returncode != 0:

            print(
                f"    attempt {attempt}: "
                f"EDirect return code "
                f"{result.returncode}"
            )

            if result.stderr.strip():

                print(
                    "    "
                    + result.stderr
                    .strip()
                    .replace(
                        "\n",
                        " "
                    )[:300]
                )

            time.sleep(
                attempt * 3
            )

            continue


        lines = [
            x
            for x in
            result.stdout.splitlines()
            if x.strip()
        ]


        has_series = False

        for line in lines:

            fields = line.split(
                "\t"
            )

            if (
                len(fields) > 0 and
                fields[0] == gse
            ):
                has_series = True
                break


        if not has_series:

            print(
                f"    attempt {attempt}: "
                "GSE record not found"
            )

            time.sleep(
                attempt * 3
            )

            continue


        tmp_file = (
            CACHE /
            f".{gse}.tmp"
        )


        with open(
            tmp_file,
            "w",
            encoding="utf-8"
        ) as f:

            f.write(
                result.stdout
            )


        tmp_file.replace(
            cache_file
        )


        success = True

        print(
            f"    success: "
            f"{len(lines)} records"
        )

        break


    if not success:

        print(
            f"    FAILED: {gse}"
        )

        failed.append(
            gse
        )


    # deliberately conservative rate
    time.sleep(
        1.0
    )


# ============================================================
# 4. Parse cached GSM records
# ============================================================

sample_rows = []
series_rows = []


for gse in targets:

    cache_file = (
        CACHE /
        f"{gse}.tsv"
    )


    if not cache_is_valid(
        cache_file,
        gse
    ):

        continue


    gsm_rows = []


    with open(
        cache_file,
        encoding="utf-8",
        errors="replace"
    ) as f:

        for line in f:

            fields = (
                line.rstrip("\n")
                .split("\t")
            )


            if not fields:
                continue


            accession = fields[0]


            if re.fullmatch(
                r"GSM[0-9]+",
                accession
            ):

                title = (
                    fields[1]
                    if len(fields) > 1
                    else ""
                )

                source = (
                    fields[2]
                    if len(fields) > 2
                    else ""
                )


                row = {
                    "GSE": gse,
                    "GSM": accession,
                    "Sample_title": title,
                    "Source_name": source
                }


                gsm_rows.append(
                    row
                )

                sample_rows.append(
                    row
                )


    if gse in series_meta.index:

        sm = series_meta.loc[
            gse
        ]

        title = sm[
            "Series_title"
        ]

        summary = sm[
            "Series_summary"
        ]

    else:

        title = ""
        summary = ""


    series_rows.append({

        "GSE": gse,

        "Series_title":
            title,

        "Series_summary":
            summary,

        "N_GSM":
            len(gsm_rows)
    })


# ============================================================
# 5. Create tables
# ============================================================

samples = pd.DataFrame(
    sample_rows,
    columns=[
        "GSE",
        "GSM",
        "Sample_title",
        "Source_name"
    ]
)


series = pd.DataFrame(
    series_rows,
    columns=[
        "GSE",
        "Series_title",
        "Series_summary",
        "N_GSM"
    ]
)


samples.to_csv(
    OUT_SAMPLES,
    sep="\t",
    index=False
)


series.to_csv(
    OUT_SERIES,
    sep="\t",
    index=False
)


# ============================================================
# 6. Collapse sample design by study
# ============================================================

def collapse_unique(
    values
):

    out = []
    seen = set()

    for value in values:

        value = str(
            value
        ).strip()

        if (
            value and
            value not in seen
        ):

            seen.add(
                value
            )

            out.append(
                value
            )

    return " || ".join(
        out
    )


def collapse_all(
    values
):

    return " || ".join(
        [
            str(x).strip()
            for x in values
            if str(x).strip()
        ]
    )


if len(samples) > 0:

    design = (
        samples
        .groupby(
            "GSE",
            as_index=False
        )
        .agg(

            N_samples=(
                "GSM",
                "size"
            ),

            Sources=(
                "Source_name",
                collapse_unique
            ),

            Sample_titles=(
                "Sample_title",
                collapse_all
            )
        )
    )

else:

    design = pd.DataFrame(
        columns=[
            "GSE",
            "N_samples",
            "Sources",
            "Sample_titles"
        ]
    )


out = series.merge(
    design,
    on="GSE",
    how="left"
)


out["N_samples"] = (
    pd.to_numeric(
        out["N_samples"],
        errors="coerce"
    )
    .fillna(0)
    .astype(int)
)


out["Sources"] = (
    out["Sources"]
    .fillna("")
)

out["Sample_titles"] = (
    out["Sample_titles"]
    .fillna("")
)


out.to_csv(
    OUT_DESIGN,
    sep="\t",
    index=False
)


# ============================================================
# 7. Failed records
# ============================================================

with open(
    OUT_FAILED,
    "w"
) as f:

    for gse in failed:

        f.write(
            gse + "\n"
        )


# ============================================================
# 8. Final report
# ============================================================

print()
print(
    "========================================"
)
print(
    "AUDIT PACK V2 COMPLETE"
)
print(
    "========================================"
)

print(
    "Requested GSE :",
    len(targets)
)

print(
    "Completed GSE :",
    len(series)
)

print(
    "Failed GSE    :",
    len(failed)
)

print(
    "Total GSM     :",
    len(samples)
)

print()
print(
    "GSE with zero direct GSM:"
)

zero = series[
    series["N_GSM"] == 0
]["GSE"].tolist()

if zero:

    for gse in zero:
        print(
            "  ",
            gse
        )

else:

    print(
        "  none"
    )


print()
print(
    "Outputs:"
)

print(
    " ",
    OUT_SERIES
)

print(
    " ",
    OUT_SAMPLES
)

print(
    " ",
    OUT_DESIGN
)

print(
    " ",
    OUT_FAILED
)
