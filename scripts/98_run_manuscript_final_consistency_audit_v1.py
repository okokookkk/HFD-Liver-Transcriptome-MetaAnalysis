#!/usr/bin/env python3
from pathlib import Path
import re
import hashlib
import pandas as pd

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()
OUTDIR = ROOT / "00_metadata/manuscript/final_consistency_audit_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT = OUTDIR / "HFD_manuscript_final_consistency_audit_v1.txt"

def rel(p):
    try:
        return str(p.relative_to(ROOT))
    except Exception:
        return str(p)

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

def read_text(path):
    try:
        return path.read_text(errors="replace")
    except Exception:
        return ""

def emit(lines, s=""):
    lines.append(str(s))

def evidence(lines, path, patterns, context=2, max_hits=40):
    if not path.exists():
        emit(lines, f"FILE_MISSING={rel(path)}")
        return 0
    txt = read_text(path).splitlines()
    regs = [re.compile(x, re.I) for x in patterns]
    hit_idx = [i for i, line in enumerate(txt) if any(r.search(line) for r in regs)]
    if not hit_idx:
        emit(lines, f"NO_MATCH={rel(path)}")
        return 0
    shown = set()
    n = 0
    for i in hit_idx:
        if n >= max_hits:
            break
        lo=max(0, i-context); hi=min(len(txt), i+context+1)
        key=(lo,hi)
        if key in shown:
            continue
        shown.add(key); n += 1
        emit(lines, f"--- {rel(path)} L{lo+1}-L{hi} ---")
        for j in range(lo,hi):
            emit(lines, f"{j+1}: {txt[j]}")
    return len(hit_idx)

lines=[]
emit(lines, "===== HFD MANUSCRIPT FINAL CONSISTENCY AUDIT v1 =====")
emit(lines)
emit(lines, f"ROOT={ROOT}")
emit(lines, "MODE=READ_ONLY")
emit(lines, "PURPOSE=MANUSCRIPT_FREEZE_NOT_NEW_INFERENCE")
emit(lines)

# ------------------------------------------------------------------
# 1. Frozen production universe
# ------------------------------------------------------------------
emit(lines, "===== 1. FROZEN PRODUCTION UNIVERSE =====")
prod = ROOT / "00_metadata/processing_production/final_audit_v1/HFD_primary64_final_production_audit_summary_v1.txt"
if prod.exists():
    emit(lines, f"FILE={rel(prod)}")
    emit(lines, f"SHA256={sha256(prod)}")
    for ln in read_text(prod).splitlines():
        if any(k in ln for k in [
            "FINAL_STATUS=", "Independent_GSE=", "Frozen_GSM=", "Frozen_SRR=",
            "Frozen_FASTQ=", "Common_gene_rows=", "Common_Geneid_SHA256=",
            "DIRECT_PRODUCTION_PASS=", "PASS_AFTER_ADAPTER_TRIM=",
            "PASS_AFTER_STRANDEDNESS_ADJUDICATION="
        ]):
            emit(lines, ln)
else:
    emit(lines, "STATUS=FAIL_MISSING_PRODUCTION_AUDIT")
emit(lines)

# ------------------------------------------------------------------
# 2. Shared GSM / sampling overlap across contrasts
# ------------------------------------------------------------------
emit(lines, "===== 2. SHARED GSM ACROSS CONTRASTS =====")
contrast_audit = ROOT / "00_metadata/analysis_inputs/frozen_v1/HFD_primary105_contrast_analysis_input_audit_v1_1.tsv"
if contrast_audit.exists():
    x = pd.read_csv(contrast_audit, sep="\t", dtype=str).fillna("")
    emit(lines, f"FILE={rel(contrast_audit)}")
    emit(lines, f"ROWS={len(x)}")
    emit(lines, f"UNIQUE_GSE={x['GSE'].nunique() if 'GSE' in x.columns else 'NA'}")
    emit(lines, f"UNIQUE_CONTRAST={x['Contrast_id'].nunique() if 'Contrast_id' in x.columns else 'NA'}")
    if "Any_shared_GSM" in x.columns:
        vc = x["Any_shared_GSM"].value_counts(dropna=False).to_dict()
        emit(lines, "Any_shared_GSM_COUNTS=" + ";".join(f"{k}:{v}" for k,v in vc.items()))
        yes = x[x["Any_shared_GSM"].str.upper().eq("YES")]
        emit(lines, f"SHARED_GSM_CONTRASTS={len(yes)}")
        if len(yes):
            cols=[c for c in ["Contrast_id","GSE","HFD_GSMs","CONTROL_GSMs",
                              "N_HFD_GSM_shared_across_contrasts",
                              "N_CONTROL_GSM_shared_across_contrasts"] if c in yes.columns]
            emit(lines, "SHARED_GSM_DETAILS:")
            emit(lines, yes[cols].to_csv(sep="\t", index=False).rstrip())
    else:
        emit(lines, "STATUS=FAIL_Any_shared_GSM_COLUMN_MISSING")
else:
    emit(lines, "STATUS=FAIL_CONTRAST_AUDIT_MISSING")
emit(lines)

# ------------------------------------------------------------------
# 3. Human frozen group sizes: primary and strict4
# ------------------------------------------------------------------
emit(lines, "===== 3. HUMAN GROUP SIZE FREEZE =====")
grp = ROOT / "00_metadata/human_validation/frozen_groups_v1_2/human_MASLD_primary_sample_groups_frozen_v1_2.tsv"
strict = ROOT / "00_metadata/human_validation/frozen_groups_v1_2/GSE130970_strict4_control_sensitivity_groups_v1_2.tsv"

def group_counts(path):
    if not path.exists():
        return None, f"MISSING={rel(path)}"
    z = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    # Try common group columns
    groupcol=None
    for c in ["Analysis_group","group","Group","Case_control","case_control","Primary_group"]:
        if c in z.columns:
            groupcol=c; break
    if groupcol is None:
        # identify a 2-level candidate with case/control-like values
        for c in z.columns:
            vals=set(v.upper() for v in z[c].dropna().astype(str).unique())
            if vals and vals <= {"CASE","CONTROL","MASLD","NORMAL","DISEASE"} and len(vals)>=2:
                groupcol=c; break
    return z, groupcol

z, gc = group_counts(grp)
if isinstance(z, pd.DataFrame):
    emit(lines, f"PRIMARY_FILE={rel(grp)}")
    emit(lines, f"PRIMARY_ROWS={len(z)}")
    if gc:
        emit(lines, f"PRIMARY_GROUP_COLUMN={gc}")
        if "GSE" in z.columns:
            for gse, zz in z.groupby("GSE"):
                emit(lines, f"{gse}=" + ";".join(f"{k}:{v}" for k,v in zz[gc].value_counts().to_dict().items()))
        else:
            emit(lines, "PRIMARY_COUNTS=" + ";".join(f"{k}:{v}" for k,v in z[gc].value_counts().to_dict().items()))
    else:
        emit(lines, "PRIMARY_GROUP_COLUMN=UNRESOLVED")
else:
    emit(lines, gc)

z2, gc2 = group_counts(strict)
if isinstance(z2, pd.DataFrame):
    emit(lines, f"STRICT4_FILE={rel(strict)}")
    emit(lines, f"STRICT4_ROWS={len(z2)}")
    if gc2:
        emit(lines, f"STRICT4_GROUP_COLUMN={gc2}")
        emit(lines, "STRICT4_COUNTS=" + ";".join(f"{k}:{v}" for k,v in z2[gc2].value_counts().to_dict().items()))
    else:
        emit(lines, "STRICT4_GROUP_COLUMN=UNRESOLVED")
else:
    emit(lines, gc2)
emit(lines)

# ------------------------------------------------------------------
# 4. Primary meta-analysis formulas / direction consistency / PI / I2
# ------------------------------------------------------------------
emit(lines, "===== 4. PRIMARY META DEFINITIONS =====")
scripts = sorted((ROOT/"scripts").glob("*"))
meta_candidates=[]
for q in scripts:
    if not q.is_file():
        continue
    txt=read_text(q)
    if ("GSE_balanced_direction_consistency" in txt or
        "prediction95_lower_new_GSE" in txt or
        "tau2_contrast_within_GSE" in txt or
        "rma.mv" in txt):
        meta_candidates.append(q)
emit(lines, f"META_SCRIPT_CANDIDATES={len(meta_candidates)}")
for q in meta_candidates:
    emit(lines, f"CANDIDATE={rel(q)}")

# Prefer newest-looking primary multilevel meta scripts.
for q in meta_candidates:
    if re.search(r"(gene|meta|multilevel)", q.name, re.I):
        emit(lines)
        emit(lines, f"### EVIDENCE FROM {rel(q)}")
        evidence(lines, q, [
            r"GSE_balanced_direction_consistency",
            r"direction_consistency",
            r"prediction95",
            r"prediction.*interval",
            r"I2_",
            r"tau2_",
            r"typical_sampling_variance",
            r"rma\.mv",
            r"random\s*=",
            r"method\s*=\s*[\"']REML"
        ], context=3, max_hits=30)
emit(lines)

# ------------------------------------------------------------------
# 5. Moderator denominators and duration+fat denominator
# ------------------------------------------------------------------
emit(lines, "===== 5. MODERATOR DENOMINATORS =====")
mod_summary_candidates=[]
for q in ROOT.rglob("*.txt"):
    sp=rel(q).lower()
    if "moderator" in sp or "duration" in sp:
        txt=read_text(q)
        if any(k in txt for k in ["2744","14588","14421","13300","2659","14195"]):
            mod_summary_candidates.append(q)
for q in sorted(mod_summary_candidates)[:30]:
    emit(lines, f"--- {rel(q)} ---")
    for ln in read_text(q).splitlines():
        if any(k in ln.lower() for k in [
            "fit", "model", "fdr", "duration", "fat", "control", "14195",
            "14421", "14588", "13300", "2744", "2659"
        ]):
            emit(lines, ln)
emit(lines)

# ------------------------------------------------------------------
# 6. Software/package versions
# ------------------------------------------------------------------
emit(lines, "===== 6. SOFTWARE / PACKAGE VERSION EVIDENCE =====")
version_patterns = [
    r"STAR[-_ ]?version", r"STAR-2\.7\.10b", r"subread-2\.0\.6",
    r"DESeq2_version", r"metafor_version", r"edgeR_version",
    r"limma_version", r"fgsea_version", r"msigdbr_version",
    r"AnnotationDbi_version", r"org\.Hs\.eg\.db.*version",
    r"R_VERSION", r"sessionInfo", r"packageVersion"
]
hits=[]
for q in list((ROOT/"scripts").glob("*")) + list((ROOT/"00_metadata").rglob("*.txt")):
    if not q.is_file():
        continue
    txt=read_text(q)
    for i,line in enumerate(txt.splitlines(),1):
        if any(re.search(p,line,re.I) for p in version_patterns):
            hits.append((rel(q),i,line.strip()))
for path,ln,line in hits[:300]:
    emit(lines, f"{path}:L{ln}: {line}")
emit(lines)

# ------------------------------------------------------------------
# 7. Hallmark source/release
# ------------------------------------------------------------------
emit(lines, "===== 7. HALLMARK SOURCE / RELEASE =====")
for q in scripts:
    if q.is_file() and re.search(r"(Hallmark|GSEA|pathway)", q.name, re.I):
        txt=read_text(q)
        if any(k.lower() in txt.lower() for k in ["msigdbr","2025.1","hallmark","fgsea"]):
            emit(lines, f"### {rel(q)}")
            evidence(lines, q, [
                r"msigdbr", r"2025\.1", r"Hallmark", r"fgsea",
                r"collection", r"species", r"category", r"subcollection"
            ], context=2, max_hits=25)
emit(lines)

# ------------------------------------------------------------------
# 8. Systematic-search provenance hunt
# ------------------------------------------------------------------
emit(lines, "===== 8. SYSTEMATIC SEARCH PROVENANCE HUNT =====")
candidate_names=[]
for q in ROOT.rglob("*"):
    if not q.is_file():
        continue
    low=q.name.lower()
    if any(k in low for k in ["prisma","systematic","search","screen","eligib"]):
        if q.suffix.lower() in {".txt",".tsv",".csv",".md",".py",".r",".sh",".json",".yaml",".yml"}:
            candidate_names.append(q)
emit(lines, f"SEARCH_FILE_CANDIDATES={len(candidate_names)}")
for q in candidate_names[:200]:
    emit(lines, rel(q))
emit(lines, "--- exact-query/date-like lines ---")
query_hits=[]
for q in candidate_names:
    txt=read_text(q)
    for i,line in enumerate(txt.splitlines(),1):
        if re.search(r"(search date|date searched|query\s*=|search query|exact query|esearch|GEO.*query|PubMed.*query)", line, re.I):
            query_hits.append((rel(q),i,line.strip()))
for path,ln,line in query_hits[:250]:
    emit(lines, f"{path}:L{ln}: {line}")
if not query_hits:
    emit(lines, "SYSTEMATIC_SEARCH_EXACT_QUERY_DATE=NOT_RECOVERED")
emit(lines)

# ------------------------------------------------------------------
# 9. Possible duplicated source cohorts across GSEs
# ------------------------------------------------------------------
emit(lines, "===== 9. POSSIBLE CROSS-GSE SOURCE-COHORT DUPLICATION SCREEN =====")
if contrast_audit.exists():
    x = pd.read_csv(contrast_audit, sep="\t", dtype=str).fillna("")
    # Gather PMID/DOI tokens from all provenance/citation/evidence columns.
    provcols=[c for c in x.columns if any(k in c.lower() for k in ["citation","source","evidence"])]
    token_rows=[]
    for _,r in x.iterrows():
        gse=r.get("GSE","")
        text=" || ".join(str(r[c]) for c in provcols)
        pmids=set(re.findall(r"PMID[:\s]*(\d{6,9})", text, re.I))
        dois=set(re.findall(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+", text, re.I))
        for t in pmids:
            token_rows.append(("PMID",t,gse))
        for t in dois:
            token_rows.append(("DOI",t.rstrip(".,;"),gse))
    if token_rows:
        d=pd.DataFrame(token_rows,columns=["type","token","GSE"]).drop_duplicates()
        g=d.groupby(["type","token"])["GSE"].agg(lambda s: sorted(set(s))).reset_index()
        dup=g[g["GSE"].map(len)>1]
        emit(lines, f"SHARED_PUBLICATION_TOKENS_ACROSS_GSE={len(dup)}")
        for _,r in dup.iterrows():
            emit(lines, f"{r['type']}={r['token']} GSE={','.join(r['GSE'])}")
        emit(lines, "NOTE=Shared publication token is a flag only; it does not prove shared animals/cohorts.")
    else:
        emit(lines, "NO_PUBLICATION_TOKENS_EXTRACTED")
else:
    emit(lines, "SKIPPED_NO_CONTRAST_AUDIT")
emit(lines)

# ------------------------------------------------------------------
# 10. Human integrated translational freeze + ortholog typo evidence
# ------------------------------------------------------------------
emit(lines, "===== 10. HUMAN TRANSLATIONAL FREEZE =====")
human = ROOT / "00_metadata/human_validation/analysis_outputs/final_summary_v1/human_validation_integrated_summary_v1.txt"
if human.exists():
    emit(lines, f"FILE={rel(human)}")
    for ln in read_text(human).splitlines():
        if any(k in ln for k in [
            "STATUS=", "Primary_common_genes=", "Final_joined_translational_genes=",
            "Overall_human_same_direction_fraction=", "Primary_OR_per_10pp=",
            "Primary_95CI=", "Primary_P=", "Mouse_vs_human_effect_Pearson_r=",
            "Hallmark_pathways_joined=", "Mouse_vs_human_RE_NES_Pearson_r="
        ]):
            emit(lines, ln)
else:
    emit(lines, "STATUS=FAIL_HUMAN_FINAL_SUMMARY_MISSING")

orthosum = ROOT / "00_metadata/human_validation/analysis_outputs/mouse_human_concordance_v1_6/mouse_human_reproducibility_concordance_summary_v1_6.txt"
if orthosum.exists():
    emit(lines, "--- ORTHOLOG PROVENANCE LABEL CHECK ---")
    for ln in read_text(orthosum).splitlines():
        if "ORTHOLOG" in ln.upper() or "ENSEMBL_RELEASE" in ln.upper():
            emit(lines, ln)
emit(lines)

# ------------------------------------------------------------------
# Final checklist
# ------------------------------------------------------------------
emit(lines, "===== FINAL MANUSCRIPT FREEZE CHECKLIST =====")
checks = [
    "SYSTEMATIC_SEARCH_EXACT_DATE_AND_QUERY",
    "SHARED_GSM_ACROSS_CONTRASTS",
    "GSE_BALANCED_DIRECTION_CONSISTENCY_EXACT_DEFINITION",
    "MULTILEVEL_I2_AND_PI_EXACT_FORMULAS",
    "STRICT4_GROUP_COUNTS_74_CASE_4_CONTROL",
    "MODERATOR_MODEL_DENOMINATORS",
    "DURATION_PLUS_FAT_MODEL_DENOMINATOR",
    "SOFTWARE_PACKAGE_VERSIONS",
    "HALLMARK_RELEASE_AND_GENESET_SOURCE",
    "POSSIBLE_SHARED_ANIMAL_COHORTS_ACROSS_GSE",
    "ORTHOLOG_PROVENANCE_LABEL_TYPO"
]
for c in checks:
    emit(lines, f"{c}=REVIEW_EVIDENCE_ABOVE")

OUT.write_text("\n".join(lines) + "\n")
print("===== HFD MANUSCRIPT FINAL CONSISTENCY AUDIT v1 =====")
print(f"Created: {OUT}")
print("MODE=READ_ONLY")
print("No analysis outputs or frozen inputs were modified.")
print()
print("Print report with:")
print(f"cat {OUT}")

