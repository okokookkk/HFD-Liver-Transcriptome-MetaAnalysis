#!/usr/bin/env python3

import csv
import io
import re
import subprocess
import time
from datetime import datetime
from pathlib import Path

import pandas as pd


BASE = Path("00_metadata/run_mapping")

INFILE = BASE / "primary107_missing_srx_gsm_v1.tsv"

OUT_SUMMARY = BASE / "primary107_missing_srx_reconciliation_v1.tsv"
OUT_RUNS = BASE / "primary107_missing_srx_rescued_runs_v1.tsv"
OUT_UNRESOLVED = BASE / "primary107_missing_srx_unresolved_v1.tsv"
OUT_LOG = BASE / "primary107_missing_srx_reconciliation_summary_v1.txt"


def extract_biosample(text):
    m = re.search(r"\b(SAM[NED]\d+)\b", str(text))
    return m.group(1) if m else ""


def fetch_runinfo(query):
    """
    Query SRA with EDirect and return RunInfo rows.
    """
    cmd = (
        f'esearch -db sra -query "{query}" '
        f'| efetch -format runinfo'
    )

    p = subprocess.run(
        ["bash", "-lc", cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=180,
    )

    if p.returncode != 0:
        return [], p.stderr.strip()

    text = p.stdout.strip()

    if not text:
        return [], ""

    try:
        reader = csv.DictReader(io.StringIO(text))
        rows = [dict(r) for r in reader]
    except Exception as e:
        return [], f"CSV parse error: {e}"

    return rows, ""


def exact_biosample_filter(rows, biosample):
    if not rows:
        return []

    if not biosample:
        return rows

    exact = []

    for r in rows:
        observed = (
            r.get("BioSample", "")
            or r.get("Sample", "")
            or ""
        ).strip()

        if observed == biosample:
            exact.append(r)

    return exact


def main():

    df = pd.read_csv(
        INFILE,
        sep="\t",
        dtype=str,
    ).fillna("")

    print("============================================")
    print("MISSING SRX RECONCILIATION V1")
    print("============================================")
    print("Input GSM:", len(df))
    print()

    reconciliation = []
    rescued_runs = []
    unresolved = []

    retrieval_time = datetime.now().astimezone().isoformat()

    for i, row in df.iterrows():

        gse = row["GSE"]
        gsm = row["GSM"]

        biosample = extract_biosample(
            row.get("All_relations", "")
        )

        print(
            f"[{i+1:02d}/{len(df)}] "
            f"{gse} {gsm} BioSample={biosample}"
        )

        matched = []
        query_method = ""
        query_value = ""
        errors = []

        # -------------------------------------------------
        # 1. Primary rescue:
        #    exact BioSample -> SRA
        # -------------------------------------------------
        if biosample:

            rows, err = fetch_runinfo(biosample)

            if err:
                errors.append(
                    f"BioSample query: {err}"
                )

            matched = exact_biosample_filter(
                rows,
                biosample,
            )

            if matched:
                query_method = "BioSample_to_SRA"
                query_value = biosample

        # -------------------------------------------------
        # 2. Fallback:
        #    GSM -> SRA, but still require exact BioSample
        #    when BioSample is available.
        # -------------------------------------------------
        if not matched:

            rows, err = fetch_runinfo(gsm)

            if err:
                errors.append(
                    f"GSM query: {err}"
                )

            matched = exact_biosample_filter(
                rows,
                biosample,
            )

            if matched:
                query_method = "GSM_to_SRA"
                query_value = gsm

        # -------------------------------------------------
        # Collect exact SRX / SRR mappings
        # -------------------------------------------------
        srx = sorted({
            r.get("Experiment", "").strip()
            for r in matched
            if r.get("Experiment", "").strip()
        })

        srr = sorted({
            r.get("Run", "").strip()
            for r in matched
            if r.get("Run", "").strip()
        })

        status = (
            "RESCUED"
            if len(srx) > 0 and len(srr) > 0
            else "UNRESOLVED"
        )

        rec = {
            "GSE": gse,
            "GSM": gsm,
            "Sample_title": row.get(
                "Sample_title", ""
            ),
            "BioSample": biosample,
            "Original_SRX": row.get("SRX", ""),
            "Resolution_status": status,
            "Resolution_method": query_method,
            "Query_accession": query_value,
            "N_SRX": len(srx),
            "SRX_list": ";".join(srx),
            "N_SRR": len(srr),
            "SRR_list": ";".join(srr),
            "Retrieval_time": retrieval_time,
            "Notes": " | ".join(errors),
        }

        reconciliation.append(rec)

        if status == "RESCUED":

            print(
                f"    RESCUED: "
                f"{len(srx)} SRX / "
                f"{len(srr)} SRR"
            )

            for r in matched:

                rescued_runs.append({
                    "GSE": gse,
                    "GSM": gsm,
                    "Sample_title":
                        row.get("Sample_title", ""),
                    "BioSample_expected": biosample,

                    "SRX":
                        r.get("Experiment", ""),
                    "SRR":
                        r.get("Run", ""),
                    "BioSample":
                        r.get("BioSample", ""),
                    "BioProject":
                        r.get("BioProject", ""),

                    "LibraryLayout":
                        r.get("LibraryLayout", ""),
                    "LibraryStrategy":
                        r.get("LibraryStrategy", ""),
                    "LibrarySource":
                        r.get("LibrarySource", ""),

                    "Platform":
                        r.get("Platform", ""),
                    "Model":
                        r.get("Model", ""),

                    "spots":
                        r.get("spots", ""),
                    "bases":
                        r.get("bases", ""),
                    "avgLength":
                        r.get("avgLength", ""),
                    "size_MB":
                        r.get("size_MB", ""),

                    "ReleaseDate":
                        r.get("ReleaseDate", ""),

                    "Resolution_method":
                        query_method,
                    "Retrieval_time":
                        retrieval_time,
                })

        else:

            print("    UNRESOLVED")

            unresolved.append(rec)

        time.sleep(0.4)

    # -----------------------------------------------------
    # Write outputs
    # -----------------------------------------------------

    pd.DataFrame(
        reconciliation
    ).to_csv(
        OUT_SUMMARY,
        sep="\t",
        index=False,
    )

    pd.DataFrame(
        rescued_runs
    ).to_csv(
        OUT_RUNS,
        sep="\t",
        index=False,
    )

    pd.DataFrame(
        unresolved
    ).to_csv(
        OUT_UNRESOLVED,
        sep="\t",
        index=False,
    )

    n_rescued = sum(
        x["Resolution_status"] == "RESCUED"
        for x in reconciliation
    )

    unique_srx = set()
    unique_srr = set()

    for x in reconciliation:

        if x["SRX_list"]:
            unique_srx.update(
                x["SRX_list"].split(";")
            )

        if x["SRR_list"]:
            unique_srr.update(
                x["SRR_list"].split(";")
            )

    with OUT_LOG.open("w") as f:

        f.write(
            "PRIMARY107 MISSING SRX "
            "RECONCILIATION V1\n"
        )
        f.write(
            "====================================\n"
        )

        f.write(
            f"Input missing-SRX GSM: "
            f"{len(df)}\n"
        )

        f.write(
            f"Rescued GSM: "
            f"{n_rescued}\n"
        )

        f.write(
            f"Unresolved GSM: "
            f"{len(unresolved)}\n"
        )

        f.write(
            f"Unique rescued SRX: "
            f"{len(unique_srx)}\n"
        )

        f.write(
            f"Unique rescued SRR: "
            f"{len(unique_srr)}\n"
        )

        f.write(
            f"Retrieval time: "
            f"{retrieval_time}\n"
        )

        f.write("\nOutputs:\n")
        f.write(f"  {OUT_SUMMARY}\n")
        f.write(f"  {OUT_RUNS}\n")
        f.write(f"  {OUT_UNRESOLVED}\n")

    print()
    print("============================================")
    print("COMPLETE")
    print("============================================")
    print("Input GSM :", len(df))
    print("Rescued   :", n_rescued)
    print("Unresolved:", len(unresolved))
    print("Unique SRX:", len(unique_srx))
    print("Unique SRR:", len(unique_srr))


if __name__ == "__main__":
    main()
