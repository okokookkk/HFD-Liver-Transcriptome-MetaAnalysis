#!/usr/bin/env python3
from pathlib import Path
import argparse
import re
import shutil
import subprocess
from datetime import datetime
import pandas as pd
import numpy as np

# ============================================================
# HFD adapter-exception downstream resolution v1
#
# Preconditions:
#   - 7 GSE were REVIEW_REQUIRED_ADAPTER in primary production.
#   - Adapter-only trimming completed with script 66 v1.1.
#   - All trimmed FASTQ passed FastQC Adapter Content.
#
# This script:
#   1) validates the trimmed FastQC gate,
#   2) aligns trimmed paired FASTQ with the SAME STAR index/options
#      used by production v1.3.2,
#   3) infers strandedness from STAR ReadsPerGene.out.tab using the
#      SAME call thresholds used in production,
#   4) stops a GSE at review if any SRR is ambiguous or a GSM contains
#      discordant SRR strand calls,
#   5) runs per-SRR featureCounts using the inferred -s value,
#   6) aggregates technical SRRs to biological GSM by exact summation,
#   7) writes canonical v1 count/QC outputs and a separate exception
#      resolution audit/status table.
#
# It DOES NOT modify:
#   - frozen manifest / eligibility / GSM/SRR/contrast membership
#   - raw FASTQ
#   - primary64 production master status v1.tsv
#
# It does NOT delete trimmed FASTQ or BAM by default.
# ============================================================

ROOT_DEFAULT = __import__("os").environ.get("HFD_PROJECT_ROOT", ".")
TARGETS = [
    'GSE145665',
    'GSE229697',
    'GSE239490',
    'GSE248297',
    'GSE285614',
    'GSE296500',
    'GSE324134',
]

ap = argparse.ArgumentParser()
ap.add_argument('--root', default=ROOT_DEFAULT)
ap.add_argument('--gse', action='append', choices=TARGETS,
                help='Run only selected GSE; may repeat. Default: all 7.')
ap.add_argument('--star-threads', type=int, default=16)
ap.add_argument('--fc-threads', type=int, default=4)
ap.add_argument('--cleanup-bam', action='store_true',
                help='After a PASS GSE, delete only its STAR BAM files. Trimmed FASTQ are retained.')
args = ap.parse_args()

ROOT = Path(args.root)
SELECTED = args.gse if args.gse else TARGETS

MANIFEST = ROOT / '00_metadata/technical_harmonization/frozen_v1/HFD_primary_ENA_R1R2_download_manifest_frozen_v1.tsv'
PRIMARY_STATUS = ROOT / '00_metadata/processing_production/HFD_primary64_production_status_v1.tsv'
TRIM_QC = ROOT / '00_metadata/processing_production/exception_review/HFD_adapter_exception_FastQC_before_after_v1_1.tsv'
TRIM_GSE_QC = ROOT / '00_metadata/processing_production/exception_review/HFD_adapter_exception_FastQC_GSE_summary_v1_1.tsv'
EXC_ROOT = ROOT / '00_metadata/processing_production/exception_review'

INDEX = ROOT / '00_metadata/reference/STAR_GRCm39_GENCODE_M39_sjdb150'
GTF = ROOT / '00_metadata/reference/GENCODE_M39_GRCm39/gencode.vM39.primary_assembly.annotation.gtf'
STAR = Path(__import__('os').environ.get('HFD_STAR', 'STAR'))
FC = Path(__import__('os').environ.get('HFD_FEATURECOUNTS', 'featureCounts'))

for q in [MANIFEST, PRIMARY_STATUS, TRIM_QC, TRIM_GSE_QC, INDEX, GTF, STAR, FC]:
    if not q.exists():
        raise RuntimeError(f'Required input/tool not found: {q}')

manifest = pd.read_csv(MANIFEST, sep='\t', dtype=str).fillna('')
status = pd.read_csv(PRIMARY_STATUS, sep='\t', dtype=str).fillna('')
trim_qc = pd.read_csv(TRIM_QC, sep='\t')
trim_gse_qc = pd.read_csv(TRIM_GSE_QC, sep='\t')

assert manifest['GSE'].nunique() == 64
assert status['GSE'].nunique() == 64

# -------------------- helpers --------------------
def now():
    return datetime.now().astimezone().isoformat()

def read_star_final(path):
    d = {}
    for line in path.read_text(errors='replace').splitlines():
        if '|' not in line:
            continue
        k, v = line.split('|', 1)
        d[k.strip()] = v.strip()
    return d

def num_from_star(d, key, integer=False):
    v = d.get(key, '').replace('%', '').strip()
    if v == '':
        return np.nan
    return int(float(v)) if integer else float(v)

def read_star_gene_counts(path):
    x = pd.read_csv(path, sep='\t', header=None,
                    names=['Geneid','unstranded','forward','reverse'])
    x = x[~x['Geneid'].astype(str).str.startswith('N_')].copy()
    for c in ['unstranded','forward','reverse']:
        x[c] = pd.to_numeric(x[c], errors='raise')
    u = int(x['unstranded'].sum())
    f = int(x['forward'].sum())
    r = int(x['reverse'].sum())
    denom = f + r
    ff = f / denom if denom else np.nan
    rf = r / denom if denom else np.nan
    return u, f, r, ff, rf

def strand_call(ff, rf):
    # Reproduces the production v1.3.2 decision regions evidenced by
    # existing STAR_strandedness_v1 tables:
    #   >=0.80 directional fraction => stranded
    #   forward fraction 0.35..0.65 => unstranded
    #   otherwise => manual review
    if ff >= 0.80:
        return 'FORWARD_STRANDED', 1
    if rf >= 0.80:
        return 'REVERSE_STRANDED', 2
    if 0.35 <= ff <= 0.65:
        return 'UNSTRANDED', 0
    return 'AMBIGUOUS_REVIEW', -1

def parse_fc_counts(path, srr):
    x = pd.read_csv(path, sep='\t', comment='#')
    meta_cols = {'Geneid','Chr','Start','End','Strand','Length'}
    count_cols = [c for c in x.columns if c not in meta_cols]
    if len(count_cols) != 1:
        raise RuntimeError(f'{srr}: featureCounts output expected 1 count column; got {count_cols}')
    out = x[['Geneid', count_cols[0]]].copy()
    out = out.rename(columns={count_cols[0]: srr})
    out[srr] = pd.to_numeric(out[srr], errors='raise').astype(np.int64)
    if (out[srr] < 0).any():
        raise RuntimeError(f'{srr}: negative featureCounts values')
    return out

def parse_fc_summary(path):
    x = pd.read_csv(path, sep='\t', index_col=0)
    if x.shape[1] != 1:
        raise RuntimeError(f'{path}: expected one summary column')
    s = pd.to_numeric(x.iloc[:,0], errors='raise')
    return s

def write_exc_status(gse, state, details=None):
    d = EXC_ROOT / gse
    d.mkdir(parents=True, exist_ok=True)
    p = d / f'{gse}_ADAPTER_EXCEPTION_DOWNSTREAM_STATUS_v1.txt'
    lines = [
        f'GSE={gse}',
        f'STATUS={state}',
        f'TIMESTAMP={now()}',
        'RESOLUTION_TYPE=ADAPTER_TRIMMED_FASTQ',
        'RAW_FASTQ_MODIFIED=NO',
        'FROZEN_TECHNICAL_ELIGIBILITY_CHANGED=NO',
        'FROZEN_GSM_SRR_CONTRAST_MEMBERSHIP_CHANGED=NO',
        'PRIMARY64_MASTER_STATUS_V1_MODIFIED=NO',
    ]
    if details:
        lines += [str(z) for z in details]
    p.write_text('\n'.join(lines) + '\n')
    return p

# -------------------- global preflight --------------------
for gse in SELECTED:
    z = status[status['GSE'].eq(gse)]
    if len(z) != 1 or z.iloc[0]['Production_status'] != 'REVIEW_REQUIRED_ADAPTER':
        raise RuntimeError(f'{gse}: primary status is not REVIEW_REQUIRED_ADAPTER')

    gq = trim_gse_qc[trim_gse_qc['GSE'].eq(gse)]
    if len(gq) != 1:
        raise RuntimeError(f'{gse}: missing/nonunique trimmed GSE FastQC summary')
    gq = gq.iloc[0]
    if str(gq['Trimmed_worst_adapter_status']) != 'pass' or int(gq['Trimmed_adapter_FAIL_n']) != 0 or int(gq['Trimmed_adapter_WARN_n']) != 0:
        raise RuntimeError(f'{gse}: trimmed adapter FastQC gate not all PASS')

resolution_rows = []

# -------------------- per-GSE processing --------------------
for gse in SELECTED:
    print('\n' + '='*100)
    print(f'{gse}: ADAPTER-TRIMMED DOWNSTREAM RESOLUTION START')
    write_exc_status(gse, 'RUNNING')

    mx = manifest[manifest['GSE'].eq(gse)].copy()
    if mx.empty:
        raise RuntimeError(f'{gse}: absent from frozen manifest')

    runmap = mx[['GSE','GSM','SRR','Contrast_id','Role']].drop_duplicates().sort_values(['GSM','SRR'])
    assert runmap.groupby('SRR')['GSM'].nunique().max() == 1
    assert runmap.groupby('GSM')['Role'].nunique().max() == 1
    assert runmap.groupby('GSM')['Contrast_id'].nunique().max() == 1
    expected_srr = sorted(runmap['SRR'].unique())
    expected_gsm = sorted(runmap['GSM'].unique())

    trim_dir = ROOT / '01_processed/adapter_trimmed_fastq' / gse
    star_dir = ROOT / '03_counts/star' / gse
    fc_dir = ROOT / '03_counts/featurecounts' / gse / 'per_srr'
    gsm_dir = ROOT / '03_counts/gsm' / gse
    starqc_dir = ROOT / '02_qc/star' / gse
    log_dir = ROOT / 'logs/production_exception_review' / gse
    prod_meta = ROOT / '00_metadata/processing_production' / gse
    exc_meta = EXC_ROOT / gse
    for d in [star_dir, fc_dir, gsm_dir, starqc_dir, log_dir, prod_meta, exc_meta]:
        d.mkdir(parents=True, exist_ok=True)

    # Exact frozen R1/R2 names, but files must come from trimmed directory.
    pair_rows = []
    for srr, q in mx.groupby('SRR', sort=True):
        comp = dict(zip(q['FASTQ_component'], q['FASTQ_filename']))
        if set(comp) != {'R1','R2'}:
            raise RuntimeError(f'{gse}/{srr}: frozen manifest is not exactly paired R1/R2')
        r1 = trim_dir / comp['R1']
        r2 = trim_dir / comp['R2']
        if not r1.exists() or not r2.exists():
            raise RuntimeError(f'{gse}/{srr}: trimmed pair missing')
        pair_rows.append((srr, r1, r2))

    # FastQC per-file gate: every expected trimmed FASTQ must be PASS.
    qz = trim_qc[trim_qc['GSE'].eq(gse)].copy()
    if len(qz) != 2 * len(expected_srr):
        raise RuntimeError(f'{gse}: trimmed FastQC row count mismatch: {len(qz)} vs {2*len(expected_srr)}')
    if not (qz['Trimmed_adapter_status'].astype(str) == 'pass').all():
        raise RuntimeError(f'{gse}: at least one trimmed FASTQ adapter module is not PASS')

    # ---------------- STAR ----------------
    mapping_rows = []
    for srr, r1, r2 in pair_rows:
        sd = star_dir / srr
        sd.mkdir(parents=True, exist_ok=True)
        bam = sd / 'Aligned.sortedByCoord.out.bam'
        final = sd / 'Log.final.out'
        gene = sd / 'ReadsPerGene.out.tab'

        if bam.exists() and final.exists() and gene.exists():
            print(f'{srr}: STAR outputs exist; validated by presence, skip')
        else:
            # Remove only incomplete STAR temp directory before safe rerun.
            tmp = sd / '_STARtmp'
            if tmp.exists():
                shutil.rmtree(tmp)
            for partial in [bam, final, gene]:
                partial.unlink(missing_ok=True)

            cmd = [
                str(STAR),
                '--runThreadN', str(args.star_threads),
                '--genomeDir', str(INDEX),
                '--readFilesIn', str(r1), str(r2),
                '--readFilesCommand', 'zcat',
                '--outFileNamePrefix', str(sd) + '/',
                '--outSAMtype', 'BAM', 'SortedByCoordinate',
                '--outSAMattributes', 'Standard',
                '--quantMode', 'GeneCounts',
            ]
            slog = log_dir / f'{srr}_STAR_adapter_trimmed_v1.log'
            print(f'{srr}: STAR')
            with slog.open('w') as fh:
                rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
            if rc != 0 or not (bam.exists() and final.exists() and gene.exists()):
                write_exc_status(gse, 'UNEXPECTED_STAR_FAILURE', [f'SRR={srr}', f'EXIT_CODE={rc}', f'LOG={slog.relative_to(ROOT)}'])
                raise RuntimeError(f'{gse}/{srr}: STAR failed; see {slog}')

        d = read_star_final(final)
        input_reads = num_from_star(d, 'Number of input reads', integer=True)
        unique_n = num_from_star(d, 'Uniquely mapped reads number', integer=True)
        unique_pct = num_from_star(d, 'Uniquely mapped reads %')
        multi_n = num_from_star(d, 'Number of reads mapped to multiple loci', integer=True)
        multi_pct = num_from_star(d, '% of reads mapped to multiple loci')
        too_many_n = num_from_star(d, 'Number of reads mapped to too many loci', integer=True)
        too_many_pct = num_from_star(d, '% of reads mapped to too many loci')
        mapping_rows.append([
            srr, input_reads, unique_n, unique_pct,
            multi_n, multi_pct, too_many_n, too_many_pct,
        ])

    mapping = pd.DataFrame(mapping_rows, columns=[
        'SRR','STAR_input_read_pairs','Uniquely_mapped_read_pairs','Uniquely_mapped_pct',
        'Multi_mapped_read_pairs','Multi_mapped_pct','Too_many_loci_read_pairs','Too_many_loci_pct'
    ]).merge(runmap[['GSM','SRR','Contrast_id','Role']], on='SRR', how='left').sort_values(['GSM','SRR'])
    mapping.to_csv(starqc_dir / f'{gse}_STAR_mapping_QC_v1.tsv', sep='\t', index=False)

    # ---------------- strandedness ----------------
    strand_rows = []
    for srr in expected_srr:
        gene = star_dir / srr / 'ReadsPerGene.out.tab'
        u, f, r, ff, rf = read_star_gene_counts(gene)
        call, sval = strand_call(ff, rf)
        strand_rows.append([srr,u,f,r,ff,rf,call,sval])

    strand = pd.DataFrame(strand_rows, columns=[
        'SRR','Gene_counts_unstranded','Gene_counts_forward','Gene_counts_reverse',
        'Forward_fraction','Reverse_fraction','Strand_call','featureCounts_s'
    ]).merge(runmap[['GSM','SRR']], on='SRR', how='left').sort_values(['GSM','SRR'])
    strand_file = starqc_dir / f'{gse}_STAR_strandedness_v1.tsv'
    strand.to_csv(strand_file, sep='\t', index=False)

    if (strand['featureCounts_s'] < 0).any():
        bad = strand[strand['featureCounts_s'] < 0]
        details = [
            f'AFFECTED_SRR={len(bad)}',
            'SRRS=' + ','.join(bad['SRR'].tolist()),
            f'STRAND_TABLE={strand_file.relative_to(ROOT)}',
        ]
        write_exc_status(gse, 'REVIEW_REQUIRED_AMBIGUOUS_STRANDEDNESS_AFTER_TRIMMING', details)
        resolution_rows.append([gse, 23, 'REVIEW_REQUIRED_AMBIGUOUS_STRANDEDNESS_AFTER_TRIMMING', len(expected_srr), len(expected_gsm)])
        print(f'{gse}: strandedness review required after trimming; downstream featureCounts not run')
        continue

    discord = strand.groupby('GSM')['featureCounts_s'].nunique()
    bad_gsm = discord[discord > 1]
    if len(bad_gsm):
        details = [
            f'AFFECTED_GSM={len(bad_gsm)}',
            'GSMS=' + ','.join(bad_gsm.index.tolist()),
            f'STRAND_TABLE={strand_file.relative_to(ROOT)}',
        ]
        write_exc_status(gse, 'REVIEW_REQUIRED_WITHIN_GSM_STRAND_DISCORDANCE_AFTER_TRIMMING', details)
        resolution_rows.append([gse, 24, 'REVIEW_REQUIRED_WITHIN_GSM_STRAND_DISCORDANCE_AFTER_TRIMMING', len(expected_srr), len(expected_gsm)])
        print(f'{gse}: within-GSM strand discordance; featureCounts not run')
        continue

    # ---------------- featureCounts per SRR ----------------
    fc_qc_rows = []
    srr_count_tables = []
    strand_by_srr = strand.set_index('SRR')['featureCounts_s'].to_dict()
    map_input = mapping.set_index('SRR')['STAR_input_read_pairs'].to_dict()

    for srr in expected_srr:
        sval = int(strand_by_srr[srr])
        bam = star_dir / srr / 'Aligned.sortedByCoord.out.bam'
        out = fc_dir / f'{srr}_featureCounts_v1.txt'
        summary = Path(str(out) + '.summary')
        flog = log_dir / f'{srr}_featureCounts_adapter_trimmed_v1.log'

        if not (out.exists() and summary.exists()):
            cmd = [
                str(FC),
                '-T', str(args.fc_threads),
                '-p', '--countReadPairs', '-B', '-C',
                '-s', str(sval),
                '-t', 'exon', '-g', 'gene_id',
                '-a', str(GTF),
                '-o', str(out),
                str(bam),
            ]
            print(f'{srr}: featureCounts -s {sval}')
            with flog.open('w') as fh:
                rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
            if rc != 0 or not (out.exists() and summary.exists()):
                write_exc_status(gse, 'UNEXPECTED_FEATURECOUNTS_FAILURE', [f'SRR={srr}', f'EXIT_CODE={rc}', f'LOG={flog.relative_to(ROOT)}'])
                raise RuntimeError(f'{gse}/{srr}: featureCounts failed; see {flog}')

        c = parse_fc_counts(out, srr)
        srr_count_tables.append(c)
        sm = parse_fc_summary(summary)
        assigned = int(sm.get('Assigned', 0))
        total = int(sm.sum())
        input_pairs = int(map_input[srr]) if not pd.isna(map_input[srr]) else np.nan
        assigned_per_input = assigned / input_pairs if input_pairs and not pd.isna(input_pairs) else np.nan
        fc_qc_rows.append([
            srr, sval, assigned, input_pairs, assigned_per_input,
            total,
            int(sm.get('Unassigned_MultiMapping', 0)),
            int(sm.get('Unassigned_NoFeatures', 0)),
            int(sm.get('Unassigned_Ambiguity', 0)),
        ])

    fcqc = pd.DataFrame(fc_qc_rows, columns=[
        'SRR','featureCounts_s','Assigned_fragments','STAR_input_read_pairs','Assigned_per_input_pair',
        'FC_summary_total_alignments','FC_Unassigned_MultiMapping_alignments','FC_Unassigned_NoFeatures','FC_Unassigned_Ambiguity'
    ]).merge(runmap[['GSM','SRR','Contrast_id','Role']], on='SRR', how='left').sort_values(['GSM','SRR'])
    fcqc.to_csv(starqc_dir / f'{gse}_featureCounts_QC_v1.tsv', sep='\t', index=False)

    # Merge SRR matrices by Geneid; exact gene identity is required.
    srr_counts = srr_count_tables[0]
    for q in srr_count_tables[1:]:
        srr_counts = srr_counts.merge(q, on='Geneid', how='inner', validate='one_to_one')
    expected_cols = ['Geneid'] + expected_srr
    if set(srr_counts.columns) != set(expected_cols):
        raise RuntimeError(f'{gse}: SRR matrix column mismatch')
    srr_counts = srr_counts[expected_cols]

    # ---------------- SRR -> GSM exact aggregation ----------------
    gsm_counts = pd.DataFrame({'Geneid': srr_counts['Geneid']})
    agg_rows = []
    for gsm, z in runmap.groupby('GSM', sort=True):
        srrs = sorted(z['SRR'].tolist())
        role = z['Role'].unique()
        contrast = z['Contrast_id'].unique()
        if len(role) != 1 or len(contrast) != 1:
            raise RuntimeError(f'{gse}/{gsm}: frozen biological mapping inconsistency')
        gsm_counts[gsm] = srr_counts[srrs].sum(axis=1).astype(np.int64)
        before = int(srr_counts[srrs].to_numpy().sum())
        after = int(gsm_counts[gsm].sum())
        agg_rows.append([
            gsm, role[0], contrast[0], len(srrs), ','.join(srrs),
            before, after, after-before, 'PASS' if after == before else 'FAIL'
        ])

    agg = pd.DataFrame(agg_rows, columns=[
        'GSM','Role','Contrast_id','SRR_count','SRRs',
        'SRR_total_assigned','GSM_total_after_sum','Difference','Aggregation_QC'
    ]).sort_values('GSM')
    agg_file = starqc_dir / f'{gse}_SRR_to_GSM_aggregation_QC_v1.tsv'
    agg.to_csv(agg_file, sep='\t', index=False)

    if not (agg['Aggregation_QC'] == 'PASS').all():
        write_exc_status(gse, 'SRR_TO_GSM_AGGREGATION_FAIL_AFTER_TRIMMING', [f'FAILED_GSM={int((agg["Aggregation_QC"] != "PASS").sum())}'])
        resolution_rows.append([gse, 25, 'SRR_TO_GSM_AGGREGATION_FAIL_AFTER_TRIMMING', len(expected_srr), len(expected_gsm)])
        continue

    # Canonical production-compatible outputs.
    final = gsm_dir / f'{gse}_GSM_counts_v1.tsv'
    sample_map = gsm_dir / f'{gse}_sample_map_v1.tsv'
    gsm_counts = gsm_counts[['Geneid'] + expected_gsm]
    gsm_counts.to_csv(final, sep='\t', index=False)
    runmap.to_csv(sample_map, sep='\t', index=False)

    # Final structural assertions.
    if list(gsm_counts.columns) != ['Geneid'] + expected_gsm:
        raise RuntimeError(f'{gse}: final GSM matrix structure failure')
    if gsm_counts['Geneid'].duplicated().any():
        raise RuntimeError(f'{gse}: duplicate Geneid in final matrix')
    if (gsm_counts.drop(columns='Geneid') < 0).any().any():
        raise RuntimeError(f'{gse}: negative final counts')

    details = [
        f'SRR={len(expected_srr)}',
        f'GSM={len(expected_gsm)}',
        f'GSM_AGGREGATION_PASS={int((agg.Aggregation_QC == "PASS").sum())}/{len(agg)}',
        f'FINAL={final.relative_to(ROOT)}',
        'STATUS=PASS_AFTER_ADAPTER_TRIM',
    ]
    exc_status_file = write_exc_status(gse, 'PASS_AFTER_ADAPTER_TRIM', details)

    # Canonical DONE marker is created only after all downstream checks pass.
    done = prod_meta / f'{gse}_PRODUCTION_DONE_v1.txt'
    done.write_text('\n'.join([
        f'GSE={gse}',
        'STATUS=PASS',
        f'TIMESTAMP={now()}',
        'RESOLVED_FROM=REVIEW_REQUIRED_ADAPTER',
        'RESOLUTION=ADAPTER_ONLY_TRIMMING_CUTADAPT_V1_1',
        'DOWNSTREAM_DRIVER=67_process_HFD_adapter_trimmed_exceptions_v1.py',
        f'FINAL={final.relative_to(ROOT)}',
        f'SAMPLE_MAP={sample_map.relative_to(ROOT)}',
        f'STAR_MAPPING_QC={(starqc_dir / f"{gse}_STAR_mapping_QC_v1.tsv").relative_to(ROOT)}',
        f'STAR_STRANDEDNESS={(starqc_dir / f"{gse}_STAR_strandedness_v1.tsv").relative_to(ROOT)}',
        f'FEATURECOUNTS_QC={(starqc_dir / f"{gse}_featureCounts_QC_v1.tsv").relative_to(ROOT)}',
        f'AGGREGATION_QC={agg_file.relative_to(ROOT)}',
        f'EXCEPTION_STATUS={exc_status_file.relative_to(ROOT)}',
        'RAW_FASTQ_MODIFIED=NO',
        'FROZEN_TECHNICAL_ELIGIBILITY_CHANGED=NO',
    ]) + '\n')

    if args.cleanup_bam:
        deleted = []
        for bam in sorted(star_dir.glob('*/Aligned.sortedByCoord.out.bam')):
            deleted.append([str(bam.relative_to(ROOT)), bam.stat().st_size, 'BAM'])
            bam.unlink()
        if deleted:
            pd.DataFrame(deleted, columns=['Deleted_path','Bytes','Type']).to_csv(
                exc_meta / f'{gse}_adapter_resolution_cleanup_manifest_v1.tsv', sep='\t', index=False
            )

    resolution_rows.append([gse, 0, 'PASS_AFTER_ADAPTER_TRIM', len(expected_srr), len(expected_gsm)])
    print(f'===== {gse} PASS AFTER ADAPTER TRIM =====')
    print('GSM matrix:', final)

# -------------------- batch status --------------------
res = pd.DataFrame(resolution_rows, columns=['GSE','Exit_code','Resolution_status','SRR','GSM'])
out = EXC_ROOT / 'HFD_adapter_exception_downstream_status_v1.tsv'
res.to_csv(out, sep='\t', index=False)

print('\n' + '='*100)
print('ADAPTER EXCEPTION DOWNSTREAM BATCH COMPLETE')
print(res.to_string(index=False))
print('\nSaved:', out.relative_to(ROOT))
print('Primary64 master status v1 was NOT modified.')

