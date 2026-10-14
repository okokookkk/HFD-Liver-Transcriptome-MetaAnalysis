#!/usr/bin/env python3

from pathlib import Path
import csv
import io
import time
import urllib.parse
import urllib.request
import pandas as pd

SRR_FILE = Path(
    "00_metadata/contrasts/"
    "HFD_primary_SRR_frozen_v1.txt"
)

OUTDIR = Path("00_metadata/ena_acquisition")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_FILES = OUTDIR / "HFD_primary_ENA_FASTQ_files_v1.tsv"
OUT_RUNS = OUTDIR / "HFD_primary_ENA_run_summary_v1.tsv"
OUT_MISSING = OUTDIR / "HFD_primary_ENA_missing_FASTQ_v1.tsv"
OUT_SUMMARY = OUTDIR / "HFD_primary_ENA_storage_summary_v1.txt"

BASE = "https://www.ebi.ac.uk/ena/portal/api/filereport"

srrs = [
    x.strip()
    for x in SRR_FILE.read_text().splitlines()
    if x.strip()
]

assert len(srrs) == 1817
assert len(set(srrs)) == 1817


def query_ena(srr, attempts=5):

    params = {
        "accession": srr,
        "result": "read_run",
        "fields":
            "run_accession,library_layout,"
            "fastq_ftp,fastq_md5,fastq_bytes",
        "format": "tsv",
    }

    url = BASE + "?" + urllib.parse.urlencode(params)

    last_error = ""

    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent":
                    "HFD-MetaAnalysis/1.0"
                }
            )

            with urllib.request.urlopen(
                req,
                timeout=60,
            ) as f:
                text = f.read().decode("utf-8")

            rows = list(
                csv.DictReader(
                    io.StringIO(text),
                    delimiter="\t",
                )
            )

            return rows, ""

        except Exception as e:
            last_error = str(e)
            time.sleep(attempt * 2)

    return [], last_error


file_rows = []
run_rows = []
missing_rows = []

for i, srr in enumerate(srrs, start=1):

    print(f"[{i:04d}/{len(srrs)}] {srr}")

    rows, err = query_ena(srr)

    if not rows:

        missing_rows.append({
            "SRR": srr,
            "Reason":
                "ENA_QUERY_FAILED_OR_NO_RECORD",
            "Error": err,
        })

        continue

    r = rows[0]

    ftp = (r.get("fastq_ftp") or "").strip()
    md5 = (r.get("fastq_md5") or "").strip()
    sizes = (r.get("fastq_bytes") or "").strip()

    layout = (
        r.get("library_layout") or ""
    ).strip()

    if not ftp:

        missing_rows.append({
            "SRR": srr,
            "Reason":
                "ENA_RECORD_HAS_NO_FASTQ",
            "Error": "",
        })

        continue

    ftp_list = ftp.split(";")
    md5_list = md5.split(";") if md5 else []
    size_list = (
        sizes.split(";")
        if sizes else []
    )

    if not (
        len(ftp_list)
        == len(md5_list)
        == len(size_list)
    ):
        missing_rows.append({
            "SRR": srr,
            "Reason":
                "ENA_FASTQ_FIELD_LENGTH_MISMATCH",
            "Error":
                f"ftp={len(ftp_list)} "
                f"md5={len(md5_list)} "
                f"bytes={len(size_list)}",
        })
        continue

    total_bytes = 0

    for j, (path, checksum, nbytes) in enumerate(
        zip(
            ftp_list,
            md5_list,
            size_list,
        ),
        start=1,
    ):

        nbytes = int(nbytes)
        total_bytes += nbytes

        url = path

        if not url.startswith(
            ("http://", "https://", "ftp://")
        ):
            url = "https://" + url

        filename = Path(
            urllib.parse.urlparse(url).path
        ).name

        file_rows.append({
            "SRR": srr,
            "LibraryLayout": layout,
            "FASTQ_index": j,
            "FASTQ_filename": filename,
            "ENA_URL": url,
            "ENA_MD5": checksum,
            "FASTQ_bytes": nbytes,
            "FASTQ_GB": nbytes / 1024**3,
        })

    run_rows.append({
        "SRR": srr,
        "LibraryLayout": layout,
        "N_FASTQ_files": len(ftp_list),
        "Total_FASTQ_bytes": total_bytes,
        "Total_FASTQ_GB":
            total_bytes / 1024**3,
    })

    # Be polite to API.
    time.sleep(0.05)


files = pd.DataFrame(file_rows)
runs = pd.DataFrame(run_rows)
missing = pd.DataFrame(missing_rows)

files.to_csv(
    OUT_FILES,
    sep="\t",
    index=False,
)

runs.to_csv(
    OUT_RUNS,
    sep="\t",
    index=False,
)

missing.to_csv(
    OUT_MISSING,
    sep="\t",
    index=False,
)

total_bytes = (
    files["FASTQ_bytes"].sum()
    if len(files)
    else 0
)

total_tb = total_bytes / 1024**4

with OUT_SUMMARY.open("w") as f:

    f.write(
        "HFD PRIMARY ENA FASTQ INVENTORY V1\n"
    )
    f.write(
        "=================================\n"
    )

    f.write(
        f"Frozen SRR requested: "
        f"{len(srrs)}\n"
    )

    f.write(
        f"SRR with ENA FASTQ: "
        f"{runs['SRR'].nunique() if len(runs) else 0}\n"
    )

    f.write(
        f"SRR missing ENA FASTQ: "
        f"{len(missing)}\n"
    )

    f.write(
        f"FASTQ.gz file count: "
        f"{len(files)}\n"
    )

    f.write(
        f"Total ENA FASTQ.gz bytes: "
        f"{total_bytes:,}\n"
    )

    f.write(
        f"Total ENA FASTQ.gz TiB: "
        f"{total_tb:.3f}\n"
    )

print()
print("==================================")
print("ENA INVENTORY COMPLETE")
print("==================================")
print("Frozen SRR:", len(srrs))
print(
    "ENA FASTQ SRR:",
    runs["SRR"].nunique()
    if len(runs) else 0
)
print("Missing SRR:", len(missing))
print("FASTQ.gz files:", len(files))
print(
    "Total FASTQ.gz:",
    f"{total_tb:.3f} TiB"
)
print()
print("Outputs:")
print(OUT_FILES)
print(OUT_RUNS)
print(OUT_MISSING)
print(OUT_SUMMARY)
