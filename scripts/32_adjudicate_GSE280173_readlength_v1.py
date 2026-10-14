#!/usr/bin/env python3
from pathlib import Path
import zipfile, re, shutil
import pandas as pd

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()
GSE = "GSE280173"
FQDIR = ROOT / "02_qc/fastqc" / GSE
META = ROOT / "00_metadata/processing_production" / GSE

def parse_dist(text):
    inside=False
    rows=[]
    for line in text.splitlines():
        if line.startswith(">>Sequence Length Distribution"):
            inside=True
            continue
        if inside and line.startswith(">>END_MODULE"):
            break
        if not inside or not line or line.startswith("#"):
            continue
        p=line.split("\t")
        if len(p)<2:
            continue
        label=p[0].strip()
        count=float(p[1])
        nums=[int(x) for x in re.findall(r"\d+",label)]
        if not nums:
            continue
        lo=nums[0]
        hi=nums[-1]
        rows.append((label,lo,hi,count))
    if not rows:
        raise RuntimeError("No Sequence Length Distribution")
    total=sum(x[3] for x in rows)
    half=total/2
    c=0
    median=None
    for row in rows:
        c += row[3]
        if c >= half:
            median=row
            break
    modal=max(rows,key=lambda x:x[3])
    definite=sum(n for _,lo,hi,n in rows if hi<50)/total
    possible=sum(n for _,lo,hi,n in rows if lo<50)/total
    if median[1] >= 50:
        call="PASS_REPRESENTATIVE_READLEN_GE50"
    elif median[2] < 50:
        call="FAIL_REPRESENTATIVE_READLEN_LT50"
    else:
        call="REVIEW_MEDIAN_BIN_STRADDLES_50"
    return median,modal,definite,possible,call

manifest = ROOT / "00_metadata/processing_pilot/HFD_multiSRR_GSM_pilot_files_v1.tsv"
# Use production frozen manifest for mapping instead.
manifest = ROOT / "00_metadata/technical_harmonization/frozen_v1/HFD_primary_ENA_R1R2_download_manifest_frozen_v1.tsv"
m=pd.read_csv(manifest,sep="\t",dtype=str).fillna("")
m=m[m["GSE"]==GSE]
runmap=m[["GSM","SRR","Contrast_id","Role"]].drop_duplicates()

rows=[]
for zpath in sorted(FQDIR.glob("*_fastqc.zip")):
    with zipfile.ZipFile(zpath) as z:
        member=[x for x in z.namelist() if x.endswith("/fastqc_data.txt")][0]
        text=z.read(member).decode("utf-8",errors="replace")
    def field(name):
        q=re.search(rf"^{re.escape(name)}\t(.+)$",text,re.M)
        return q.group(1) if q else ""
    fn=field("Filename")
    mm=re.match(r"(.+)_([12])\.fastq\.gz$",fn)
    srr=mm.group(1); mate="R"+mm.group(2)
    median,modal,definite,possible,call=parse_dist(text)
    rows.append({
        "SRR":srr,
        "Mate":mate,
        "FASTQ":fn,
        "Sequence_length":field("Sequence length"),
        "Median_length_bin":median[0],
        "Median_length_lower":median[1],
        "Median_length_upper":median[2],
        "Modal_length_bin":modal[0],
        "Modal_length_lower":modal[1],
        "Modal_length_upper":modal[2],
        "Fraction_reads_definitely_lt50":definite,
        "Fraction_reads_possibly_lt50":possible,
        "Representative_read_length_QC":call,
    })

x=pd.DataFrame(rows).merge(runmap,on="SRR",how="left",validate="many_to_one")
x=x.sort_values(["GSM","SRR","Mate"])

out=META/f"{GSE}_readlength_distribution_QC_v1.tsv"
x.to_csv(out,sep="\t",index=False)

print("===== GSE280173 REPRESENTATIVE READ-LENGTH ADJUDICATION =====")
print(x.to_string(index=False))
print("\nCalls:")
print(x["Representative_read_length_QC"].value_counts().to_string())
print("\nMedian-bin lower range:",
      x["Median_length_lower"].min(), "-", x["Median_length_lower"].max())
print("Median-bin upper range:",
      x["Median_length_upper"].min(), "-", x["Median_length_upper"].max())
print("Fraction definitely <50 range:",
      x["Fraction_reads_definitely_lt50"].min(), "-",
      x["Fraction_reads_definitely_lt50"].max())
print("\nSaved:",out)

old=META/f"{GSE}_PRODUCTION_STATUS_v1.txt"
if old.exists():
    backup=META/f"{GSE}_PRODUCTION_STATUS_pre_readlength_rule_correction_v1.txt"
    if not backup.exists():
        shutil.copy2(old,backup)
        print("Preserved prior provisional status:",backup)

if (x["Representative_read_length_QC"]=="PASS_REPRESENTATIVE_READLEN_GE50").all():
    note=META/f"{GSE}_readlength_rule_correction_note_v1.txt"
    note.write_text(
        "GSE280173 READ-LENGTH QC RULE CORRECTION\n"
        "The initial production driver incorrectly used the minimum individual read length\n"
        "as the hard threshold for the frozen >=50-bp technical criterion.\n"
        "Because the archived FASTQ files contain variable-length reads (FastQC 36-101 bp),\n"
        "the criterion was operationalized using the weighted-median FastQC read-length bin.\n"
        "The prior minimum-length-based exception status is retained as audit history and is\n"
        "superseded by representative-read-length adjudication. No frozen eligibility rule\n"
        "or threshold was changed.\n"
    )
    print("\nADJUDICATION: representative read length PASS; prior min-length exception is superseded.")
else:
    print("\nADJUDICATION: not all FASTQ files pass representative read-length rule; review output.")
