#!/usr/bin/env python3
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse, hashlib, gzip, re, shutil, subprocess, sys, zipfile, time
from datetime import datetime
import pandas as pd
import numpy as np

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()
MANIFEST = ROOT / '00_metadata/technical_harmonization/frozen_v1/HFD_primary_ENA_R1R2_download_manifest_frozen_v1.tsv'
INDEX = ROOT / '00_metadata/reference/STAR_GRCm39_GENCODE_M39_sjdb150'
GTF = ROOT / '00_metadata/reference/GENCODE_M39_GRCm39/gencode.vM39.primary_assembly.annotation.gtf'
STAR = Path(__import__('os').environ.get('HFD_STAR', 'STAR'))
FC = Path(__import__('os').environ.get('HFD_FEATURECOUNTS', 'featureCounts'))

ap = argparse.ArgumentParser()
ap.add_argument('--gse', required=True)
ap.add_argument('--download-workers', type=int, default=4)
ap.add_argument('--star-threads', type=int, default=16)
ap.add_argument('--fastqc-threads', type=int, default=16)
ap.add_argument('--fc-threads', type=int, default=4)
ap.add_argument('--cleanup', action='store_true')
args = ap.parse_args()
GSE = args.gse

os_paths = {
    'raw': ROOT / '01_raw/fastq' / GSE,
    'fastqc': ROOT / '02_qc/fastqc' / GSE,
    'multiqc': ROOT / '02_qc/multiqc' / GSE,
    'starqc': ROOT / '02_qc/star' / GSE,
    'star': ROOT / '03_counts/star' / GSE,
    'fc': ROOT / '03_counts/featurecounts' / GSE / 'per_srr',
    'gsm': ROOT / '03_counts/gsm' / GSE,
    'meta': ROOT / '00_metadata/processing_production' / GSE,
    'log': ROOT / 'logs/production' / GSE,
}
for p in os_paths.values(): p.mkdir(parents=True, exist_ok=True)

DONE = os_paths['meta'] / f'{GSE}_PRODUCTION_DONE_v1.txt'
FINAL = os_paths['gsm'] / f'{GSE}_GSM_counts_v1.tsv'
STATUS_FILE = os_paths['meta'] / f'{GSE}_PRODUCTION_STATUS_v1.txt'

def write_status(status, details=None):
    lines = [
        f'GSE={GSE}',
        f'STATUS={status}',
        f'TIMESTAMP={datetime.now().astimezone().isoformat()}',
    ]
    if details:
        lines.extend(str(x) for x in details)
    STATUS_FILE.write_text('\n'.join(lines) + '\n')

def processing_stop(status, code, details=None):
    write_status(status, details)
    exc = os_paths['meta'] / f'{GSE}_{status}_v1.txt'
    exc.write_text(STATUS_FILE.read_text())
    print(f'{GSE}: {status}')
    raise SystemExit(code)

write_status('RUNNING')
if DONE.exists() and FINAL.exists():
    if args.cleanup:
        deleted=[]
        for p in sorted(os_paths['raw'].glob('*.fastq.gz')):
            deleted.append([str(p.relative_to(ROOT)),p.stat().st_size,'FASTQ']); p.unlink()
        for p in sorted(os_paths['star'].glob('*/Aligned.sortedByCoord.out.bam')):
            deleted.append([str(p.relative_to(ROOT)),p.stat().st_size,'BAM']); p.unlink()
        if deleted:
            pd.DataFrame(deleted,columns=['Deleted_path','Bytes','Type']).to_csv(os_paths['meta']/f'{GSE}_cleanup_manifest_v1.tsv',sep='\t',index=False)
    print(f'{GSE}: already complete; verified final output exists')
    sys.exit(0)

x = pd.read_csv(MANIFEST, sep='\t', dtype=str).fillna('')
x = x[x['GSE'] == GSE].copy()
if x.empty: raise SystemExit(f'{GSE}: not found in frozen technical manifest')
x['FASTQ_bytes'] = pd.to_numeric(x['FASTQ_bytes'], errors='raise')
assert not x['FASTQ_filename'].duplicated().any()
runmap = x[['GSE','GSM','SRR','Contrast_id','Role']].drop_duplicates()
assert runmap['SRR'].nunique() == x['SRR'].nunique()
assert runmap.groupby('SRR')['GSM'].nunique().max() == 1
assert runmap.groupby('GSM')['Role'].nunique().max() == 1
assert runmap.groupby('GSM')['Contrast_id'].nunique().max() == 1
pair_components = x.groupby('SRR')['FASTQ_component'].apply(lambda z: set(z))
assert all(v == {'R1','R2'} for v in pair_components)
expected_srr = sorted(runmap['SRR'].unique())
expected_gsm = sorted(runmap['GSM'].unique())

print(f'===== {GSE} PRODUCTION START =====')
print('GSM:', len(expected_gsm), 'SRR:', len(expected_srr), 'FASTQ:', len(x))

# Protocol provenance from local GEO SOFT, if present
soft = ROOT / '00_metadata/geo_soft' / f'{GSE}_family.soft.gz'
if soft.exists():
    pattern = r'^\!Series_(title|summary|overall_design)|^\!Sample_(title|source_name_ch1|organism_ch1|characteristics_ch1|extract_protocol_ch1|library_strategy|library_source|library_selection|data_processing)'
    r = subprocess.run(['zgrep','-E',pattern,str(soft)], text=True, capture_output=True)
    lines = sorted(set([z for z in r.stdout.splitlines() if z]))
    (os_paths['meta'] / f'{GSE}_protocol_review_v1.txt').write_text('\n'.join(lines)+'\n')

# -------------------- acquisition --------------------
def md5sum(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(16*1024*1024), b''): h.update(b)
    return h.hexdigest()

def gzip_test(path):
    with gzip.open(path, 'rb') as f:
        while f.read(16*1024*1024): pass

def urlfix(url):
    return url if url.startswith(('http://','https://')) else 'https://' + url.lstrip('/')

def acquire(row):
    fn = row['FASTQ_filename']
    dst = os_paths['raw'] / fn
    es = int(row['FASTQ_bytes'])
    em = row['ENA_MD5'].lower()
    https_url = urlfix(row['ENA_URL'])

    # ENA fastq_ftp-style URLs point to ftp.sra.ebi.ac.uk.  Keep an FTP
    # transport fallback because transient HTTPS/proxy failures can occur.
    raw_url = row['ENA_URL'].strip()
    if raw_url.startswith('ftp://'):
        ftp_url = raw_url
    elif raw_url.startswith(('http://', 'https://')):
        ftp_url = re.sub(r'^https?://', 'ftp://', raw_url)
    else:
        ftp_url = 'ftp://' + raw_url.lstrip('/')

    attempt_rows = []
    attempt_file = os_paths['meta'] / f'{fn}.download_attempts_v1.tsv'

    def log_attempt(method, attempt, rc, before, after, note=''):
        attempt_rows.append({
            'Timestamp': datetime.now().astimezone().isoformat(),
            'FASTQ_filename': fn,
            'Method': method,
            'Attempt': attempt,
            'Exit_code': rc,
            'Bytes_before': before,
            'Bytes_after': after,
            'Expected_bytes': es,
            'Bytes_added': after - before,
            'Note': note,
        })
        pd.DataFrame(attempt_rows).to_csv(
            attempt_file, sep='\t', index=False
        )

    def validate_if_complete():
        """
        Return (valid, observed_bytes, observed_md5, gzip_qc).
        Only calculate MD5/gzip for a file at the exact expected size.
        """
        obs = dst.stat().st_size if dst.exists() else 0
        if obs != es:
            return False, obs, '', 'NOT_RUN'

        am = md5sum(dst)
        if am != em:
            return False, obs, am, 'NOT_RUN'

        try:
            gzip_test(dst)
        except Exception:
            return False, obs, am, 'FAIL'

        return True, obs, am, 'PASS'

    # Existing exact-size file is accepted only after archive MD5 + gzip.
    if dst.exists():
        if dst.stat().st_size > es:
            dst.unlink()
        elif dst.stat().st_size == es:
            valid, obs, am, gz = validate_if_complete()
            if valid:
                return [
                    row['GSM'], row['SRR'], row['FASTQ_component'], fn,
                    es, obs, em, am,
                    'EXISTING_VALID', 'PASS', 'PASS', 'PASS'
                ]
            # Exact-size but checksum/gzip-invalid files must never be resumed.
            dst.unlink(missing_ok=True)

    def curl_attempt(url, method, extra_args, outer_attempts):
        previous_size = dst.stat().st_size if dst.exists() else 0
        last_rc = None

        for attempt in range(1, outer_attempts + 1):
            if dst.exists() and dst.stat().st_size > es:
                dst.unlink()
                previous_size = 0

            before = dst.stat().st_size if dst.exists() else 0

            cmd = [
                'curl', '--silent', '--show-error', '-L', '--fail',
                '--retry', '5', '--retry-delay', '5',
                '--connect-timeout', '30',
                '--continue-at', '-',
            ] + list(extra_args) + [
                '--output', str(dst), url
            ]

            rc = subprocess.run(cmd).returncode
            last_rc = rc
            after = dst.stat().st_size if dst.exists() else 0
            log_attempt(method, attempt, rc, before, after)

            valid, obs, am, gz = validate_if_complete()
            if valid:
                return True, rc, obs, am, method

            # Exact expected size but bad checksum/gzip: restart fresh.
            if obs == es:
                log_attempt(
                    method, attempt, rc, obs, obs,
                    note=f'FULL_SIZE_VALIDATION_FAILED md5={am or "NA"} gzip={gz}'
                )
                dst.unlink(missing_ok=True)
                previous_size = 0
            else:
                grew = after > previous_size
                previous_size = after
                if attempt < outer_attempts:
                    time.sleep(min(5 * attempt, 30) if grew else 30)

        obs = dst.stat().st_size if dst.exists() else 0
        return False, last_rc, obs, '', method

    def wget_attempt(url, outer_attempts=6):
        if shutil.which('wget') is None:
            return False, None, dst.stat().st_size if dst.exists() else 0, '', 'WGET_HTTPS_UNAVAILABLE'

        last_rc = None
        for attempt in range(1, outer_attempts + 1):
            if dst.exists() and dst.stat().st_size > es:
                dst.unlink()

            before = dst.stat().st_size if dst.exists() else 0
            cmd = [
                'wget',
                '--continue',
                '--tries=3',
                '--timeout=30',
                '--read-timeout=60',
                '--output-document', str(dst),
                url
            ]

            rc = subprocess.run(cmd).returncode
            last_rc = rc
            after = dst.stat().st_size if dst.exists() else 0
            log_attempt('WGET_HTTPS', attempt, rc, before, after)

            valid, obs, am, gz = validate_if_complete()
            if valid:
                return True, rc, obs, am, 'WGET_HTTPS'

            if obs == es:
                log_attempt(
                    'WGET_HTTPS', attempt, rc, obs, obs,
                    note=f'FULL_SIZE_VALIDATION_FAILED md5={am or "NA"} gzip={gz}'
                )
                dst.unlink(missing_ok=True)

            if attempt < outer_attempts:
                time.sleep(min(10 * attempt, 30))

        obs = dst.stat().st_size if dst.exists() else 0
        return False, last_rc, obs, '', 'WGET_HTTPS'

    # Transport sequence:
    # 1) regular ENA HTTPS
    # 2) HTTPS forced to HTTP/1.1 (useful for some proxy/server disconnects)
    # 3) wget HTTPS (different client/resume implementation)
    # 4) ENA FTP transport
    strategies = [
        ('CURL_HTTPS', https_url, [], 8),
        ('CURL_HTTPS_HTTP1_1', https_url, ['--http1.1'], 8),
    ]

    for method, url, extra, ntry in strategies:
        ok, rc, obs, am, used = curl_attempt(url, method, extra, ntry)
        if ok:
            return [
                row['GSM'], row['SRR'], row['FASTQ_component'], fn,
                es, obs, em, am,
                f'DOWNLOADED_AFTER_{used}', 'PASS', 'PASS', 'PASS'
            ]

    ok, rc, obs, am, used = wget_attempt(https_url, outer_attempts=6)
    if ok:
        return [
            row['GSM'], row['SRR'], row['FASTQ_component'], fn,
            es, obs, em, am,
            f'DOWNLOADED_AFTER_{used}', 'PASS', 'PASS', 'PASS'
        ]

    ok, rc, obs, am, used = curl_attempt(
        ftp_url, 'CURL_FTP', ['--ftp-pasv'], 8
    )
    if ok:
        return [
            row['GSM'], row['SRR'], row['FASTQ_component'], fn,
            es, obs, em, am,
            f'DOWNLOADED_AFTER_{used}', 'PASS', 'PASS', 'PASS'
        ]

    # Exhausted all robust acquisition routes.  Do NOT call this a biological
    # or technical exclusion.  Preserve the partial file for future resume and
    # return a structured failure row so the GSE can be deferred.
    obs = dst.stat().st_size if dst.exists() else 0
    am = md5sum(dst) if obs == es else ''
    return [
        row['GSM'], row['SRR'], row['FASTQ_component'], fn,
        es, obs, em, am,
        'DEFERRED_ACQUISITION_FAILURE',
        'PASS' if obs == es else 'FAIL',
        'PASS' if am == em and am else ('FAIL' if am else 'NOT_RUN'),
        'NOT_RUN'
    ]

acq_rows=[]
worker_errors=[]
with ThreadPoolExecutor(max_workers=args.download_workers) as pool:
    futs={pool.submit(acquire,row): row['FASTQ_filename'] for _,row in x.iterrows()}
    for fut in as_completed(futs):
        fn = futs[fut]
        try:
            acq_rows.append(fut.result())
        except Exception as e:
            worker_errors.append((fn, repr(e)))

if worker_errors:
    pd.DataFrame(
        worker_errors,
        columns=['FASTQ_filename','Worker_exception']
    ).to_csv(
        os_paths['meta']/f'{GSE}_acquisition_worker_exceptions_v1.tsv',
        sep='\t',
        index=False
    )
    processing_stop(
        'DEFERRED_ACQUISITION_FAILURE',
        27,
        [
            f'WORKER_EXCEPTIONS={len(worker_errors)}',
            'FROZEN_TECHNICAL_ELIGIBILITY_CHANGED=NO',
            'FINAL_HARMONIZED_RAW_INCLUDED=PENDING_RETRY',
            'INTERPRETATION=ARCHIVE_TRANSFER_DEFERRED_NOT_EXCLUDED'
        ]
    )

acq = pd.DataFrame(
    acq_rows,
    columns=[
        'GSM','SRR','FASTQ_component','FASTQ_filename',
        'Expected_bytes','Observed_bytes','Expected_MD5','Observed_MD5',
        'Download_status','Size_QC','MD5_QC','Gzip_QC'
    ]
).sort_values(['GSM','SRR','FASTQ_component'])

assert len(acq) == len(x)

acq.to_csv(
    os_paths['meta']/f'{GSE}_FASTQ_acquisition_QC_v1.tsv',
    sep='\t',
    index=False
)

failed_acq = acq[
    (acq['Size_QC'] != 'PASS')
    | (acq['MD5_QC'] != 'PASS')
    | (acq['Gzip_QC'] != 'PASS')
].copy()

if len(failed_acq):
    failed_acq.to_csv(
        os_paths['meta']/f'{GSE}_acquisition_deferred_files_v1.tsv',
        sep='\t',
        index=False
    )
    processing_stop(
        'DEFERRED_ACQUISITION_FAILURE',
        27,
        [
            f'AFFECTED_FASTQ={len(failed_acq)}',
            f'AFFECTED_SRR={failed_acq["SRR"].nunique()}',
            'PARTIAL_FILES_PRESERVED_FOR_RESUME=YES',
            'FROZEN_TECHNICAL_ELIGIBILITY_CHANGED=NO',
            'FINAL_HARMONIZED_RAW_INCLUDED=PENDING_RETRY',
            'INTERPRETATION=ARCHIVE_TRANSFER_DEFERRED_NOT_EXCLUDED'
        ]
    )

# -------------------- read-length distribution helpers --------------------
def parse_fastqc_length_distribution(text):
    """Parse FastQC Sequence Length Distribution and return weighted-median metrics."""
    inside = False
    rows = []
    for line in text.splitlines():
        if line.startswith('>>Sequence Length Distribution'):
            inside = True
            continue
        if inside and line.startswith('>>END_MODULE'):
            break
        if not inside or not line or line.startswith('#'):
            continue
        p = line.split('\t')
        if len(p) < 2:
            continue
        label = p[0].strip()
        count = float(p[1])
        nums = [int(x) for x in re.findall(r'\d+', label)]
        if not nums:
            continue
        lo = nums[0]
        hi = nums[-1]
        rows.append((label, lo, hi, count))

    if not rows:
        raise RuntimeError('FastQC Sequence Length Distribution module not found')

    total = sum(r[3] for r in rows)
    if total <= 0:
        raise RuntimeError('FastQC Sequence Length Distribution has zero total count')

    half = total / 2.0
    cumulative = 0.0
    median_row = None
    for row in rows:
        cumulative += row[3]
        if cumulative >= half:
            median_row = row
            break

    modal_row = max(rows, key=lambda r: r[3])
    definitely_lt50 = sum(count for _, lo, hi, count in rows if hi < 50)
    possibly_lt50 = sum(count for _, lo, hi, count in rows if lo < 50)

    return {
        'Length_distribution_total': total,
        'Median_length_bin': median_row[0],
        'Median_length_lower': median_row[1],
        'Median_length_upper': median_row[2],
        'Modal_length_bin': modal_row[0],
        'Modal_length_lower': modal_row[1],
        'Modal_length_upper': modal_row[2],
        'Fraction_reads_definitely_lt50': definitely_lt50 / total,
        'Fraction_reads_possibly_lt50': possibly_lt50 / total,
    }


def representative_read_length_call(lower, upper):
    """Apply frozen >=50-bp criterion to the weighted-median FastQC length bin."""
    if lower >= 50:
        return 'PASS_REPRESENTATIVE_READLEN_GE50'
    if upper < 50:
        return 'FAIL_REPRESENTATIVE_READLEN_LT50'
    return 'REVIEW_MEDIAN_BIN_STRADDLES_50'

# -------------------- FastQC --------------------
fastq_files = [str(os_paths['raw']/fn) for fn in x['FASTQ_filename']]
expected_zips = [os_paths['fastqc']/(Path(fn).name.replace('.fastq.gz','_fastqc.zip')) for fn in x['FASTQ_filename']]
if not all(p.exists() for p in expected_zips):
    cmd=['fastqc','-t',str(args.fastqc_threads),'-o',str(os_paths['fastqc'])] + fastq_files
    with open(os_paths['log']/f'{GSE}_FastQC_raw_v1.log','w') as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
subprocess.run(['multiqc',str(os_paths['fastqc']),'-o',str(os_paths['multiqc']),'-n',f'{GSE}_raw_MultiQC_v1','--force'], stdout=open(os_paths['log']/f'{GSE}_MultiQC_raw_v1.log','w'), stderr=subprocess.STDOUT, check=True)

fq_rows=[]
for zpath in sorted(os_paths['fastqc'].glob('*_fastqc.zip')):
    with zipfile.ZipFile(zpath) as z:
        data_member=[n for n in z.namelist() if n.endswith('/fastqc_data.txt')][0]
        summary_member=[n for n in z.namelist() if n.endswith('/summary.txt')][0]
        text=z.read(data_member).decode('utf-8', errors='replace')
        summary=z.read(summary_member).decode('utf-8', errors='replace')
    def field(name):
        m=re.search(rf'^{re.escape(name)}\t(.+)$', text, re.M); return m.group(1) if m else ''
    fn=field('Filename'); mm=re.match(r'(.+)_([12])\.fastq\.gz$',fn)
    if not mm: continue
    srr=mm.group(1); mate='R'+mm.group(2); nseq=int(field('Total Sequences').replace(',','')); slen=field('Sequence length')
    nums=[int(v) for v in re.findall(r'\d+',slen)]
    modules={}
    for line in summary.strip().splitlines():
        p=line.split('\t')
        if len(p)>=2: modules[p[1]]=p[0]

    ld = parse_fastqc_length_distribution(text)
    rep_call = representative_read_length_call(
        ld['Median_length_lower'],
        ld['Median_length_upper']
    )

    fq_rows.append([
        srr,mate,fn,nseq,slen,min(nums),max(nums),field('%GC'),
        ld['Median_length_bin'],ld['Median_length_lower'],ld['Median_length_upper'],
        ld['Modal_length_bin'],ld['Modal_length_lower'],ld['Modal_length_upper'],
        ld['Fraction_reads_definitely_lt50'],ld['Fraction_reads_possibly_lt50'],
        rep_call,
        modules.get('Per base sequence quality',''),
        modules.get('Per sequence quality scores',''),
        modules.get('Adapter Content',''),
        modules.get('Overrepresented sequences',''),
        modules.get('Sequence Duplication Levels','')
    ])
fq=pd.DataFrame(
    fq_rows,
    columns=[
        'SRR','Mate','FASTQ','Total_sequences','Sequence_length',
        'Min_read_length','Max_read_length','Percent_GC',
        'Median_length_bin','Median_length_lower','Median_length_upper',
        'Modal_length_bin','Modal_length_lower','Modal_length_upper',
        'Fraction_reads_definitely_lt50','Fraction_reads_possibly_lt50',
        'Representative_read_length_QC',
        'Per_base_sequence_quality','Per_sequence_quality_scores',
        'Adapter_Content','Overrepresented_sequences','Sequence_Duplication_Levels'
    ]
)
fq=fq.merge(runmap[['GSM','SRR','Contrast_id','Role']], on='SRR', how='left', validate='many_to_one').sort_values(['GSM','SRR','Mate'])
assert len(fq)==2*len(expected_srr) and fq['SRR'].nunique()==len(expected_srr)
pivot=fq.pivot(index='SRR',columns='Mate',values='Total_sequences'); assert {'R1','R2'} <= set(pivot.columns)
pivot['Difference']=pivot['R1']-pivot['R2']; pivot['Pair_count_equal']=pivot['Difference'].eq(0)

# Save structural provenance BEFORE any QC stop.
fq.to_csv(os_paths['meta']/f'{GSE}_FastQC_basic_statistics_v1.tsv', sep='\t', index=False)
pivot.reset_index().to_csv(os_paths['meta']/f'{GSE}_R1R2_pair_count_QC_v1.tsv', sep='\t', index=False)

if not pivot['Pair_count_equal'].all():
    bad = pivot.loc[~pivot['Pair_count_equal']].reset_index()
    bad.to_csv(os_paths['meta']/f'{GSE}_R1R2_count_mismatch_affected_v1.tsv', sep='\t', index=False)
    processing_stop(
        'PROCESSING_QC_EXCEPTION_R1R2_COUNT_MISMATCH',
        21,
        [
            f'AFFECTED_SRR={len(bad)}',
            'FROZEN_TECHNICAL_ELIGIBILITY_CHANGED=NO',
            'FINAL_HARMONIZED_RAW_INCLUDED=NO_PENDING_ADJUDICATION'
        ]
    )

# Operationalize the frozen >=50-bp criterion using the weighted-median
# FastQC read-length distribution, NOT the minimum individual read length.
# Variable-length libraries can legitimately contain a minority of reads <50 bp.

readlen_qc = fq[
    [
        'GSM','SRR','Mate','FASTQ','Sequence_length',
        'Min_read_length','Max_read_length',
        'Median_length_bin','Median_length_lower','Median_length_upper',
        'Modal_length_bin','Modal_length_lower','Modal_length_upper',
        'Fraction_reads_definitely_lt50','Fraction_reads_possibly_lt50',
        'Representative_read_length_QC'
    ]
].copy()

readlen_qc.to_csv(
    os_paths['meta']/f'{GSE}_readlength_distribution_QC_v1.tsv',
    sep='\t',
    index=False
)

fail = fq[
    fq['Representative_read_length_QC']
    == 'FAIL_REPRESENTATIVE_READLEN_LT50'
].copy()

review = fq[
    fq['Representative_read_length_QC']
    == 'REVIEW_MEDIAN_BIN_STRADDLES_50'
].copy()

if len(fail):
    fail[
        [
            'GSM','SRR','Mate','FASTQ','Sequence_length',
            'Median_length_bin','Median_length_lower','Median_length_upper',
            'Fraction_reads_definitely_lt50','Fraction_reads_possibly_lt50'
        ]
    ].to_csv(
        os_paths['meta']/f'{GSE}_representative_readlength_LT50_affected_v1.tsv',
        sep='\t',
        index=False
    )
    processing_stop(
        'PROCESSING_QC_EXCEPTION_REPRESENTATIVE_READLEN_LT50',
        20,
        [
            f'AFFECTED_FASTQ={len(fail)}',
            f'AFFECTED_SRR={fail["SRR"].nunique()}',
            'READ_LENGTH_RULE=WEIGHTED_MEDIAN_FASTQC_LENGTH_BIN',
            'THRESHOLD_BP=50',
            'FROZEN_TECHNICAL_ELIGIBILITY_CHANGED=NO',
            'FINAL_HARMONIZED_RAW_INCLUDED=NO',
            'INTERPRETATION=POST_FREEZE_RAW_PROCESSING_EXCEPTION'
        ]
    )

if len(review):
    review[
        [
            'GSM','SRR','Mate','FASTQ','Sequence_length',
            'Median_length_bin','Median_length_lower','Median_length_upper',
            'Fraction_reads_definitely_lt50','Fraction_reads_possibly_lt50'
        ]
    ].to_csv(
        os_paths['meta']/f'{GSE}_readlength_borderline_review_v1.tsv',
        sep='\t',
        index=False
    )
    processing_stop(
        'REVIEW_REQUIRED_READ_LENGTH_BORDERLINE',
        26,
        [
            f'AFFECTED_FASTQ={len(review)}',
            f'AFFECTED_SRR={review["SRR"].nunique()}',
            'READ_LENGTH_RULE=WEIGHTED_MEDIAN_FASTQC_LENGTH_BIN',
            'THRESHOLD_BP=50',
            'FROZEN_TECHNICAL_ELIGIBILITY_CHANGED=NO'
        ]
    )

if (fq['Adapter_Content']=='FAIL').any():
    bad = fq.loc[
        fq['Adapter_Content'] == 'FAIL',
        ['GSM','SRR','Mate','FASTQ','Adapter_Content']
    ].copy()
    bad.to_csv(os_paths['meta']/f'{GSE}_adapter_FAIL_affected_v1.tsv', sep='\t', index=False)
    processing_stop(
        'REVIEW_REQUIRED_ADAPTER',
        22,
        [
            f'AFFECTED_FASTQ={len(bad)}',
            'NO_AUTOMATIC_TRIMMING_PERFORMED=YES',
            'FROZEN_TECHNICAL_ELIGIBILITY_CHANGED=NO'
        ]
    )

# -------------------- STAR --------------------
for srr in expected_srr:
    r1=os_paths['raw']/f'{srr}_1.fastq.gz'; r2=os_paths['raw']/f'{srr}_2.fastq.gz'; out=os_paths['star']/srr; out.mkdir(parents=True, exist_ok=True)
    bam=out/'Aligned.sortedByCoord.out.bam'; final=out/'Log.final.out'; gene=out/'ReadsPerGene.out.tab'
    if bam.exists() and final.exists() and gene.exists(): continue
    cmd=[str(STAR),'--runThreadN',str(args.star_threads),'--genomeDir',str(INDEX),'--readFilesIn',str(r1),str(r2),'--readFilesCommand','zcat','--outFileNamePrefix',str(out)+'/', '--outSAMtype','BAM','SortedByCoordinate','--outSAMattributes','Standard','--quantMode','GeneCounts']
    with open(os_paths['log']/f'{srr}_STAR_v1.log','w') as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)

# -------------------- STAR QC + strandedness --------------------
map_rows=[]; strand_rows=[]
for srr in expected_srr:
    log=os_paths['star']/srr/'Log.final.out'; d={}
    for line in log.read_text().splitlines():
        if '|' in line:
            k,v=line.split('|',1); d[k.strip()]=v.strip()
    map_rows.append([srr,int(d['Number of input reads']),int(d['Average input read length']),int(d['Uniquely mapped reads number']),float(d['Uniquely mapped reads %'].rstrip('%')),float(d['Average mapped length']),float(d['Mismatch rate per base, %'].rstrip('%')),float(d['% of reads mapped to multiple loci'].rstrip('%')),float(d['% of reads mapped to too many loci'].rstrip('%')),float(d['% of reads unmapped: too short'].rstrip('%')),float(d['% of reads unmapped: other'].rstrip('%'))])
    g=pd.read_csv(os_paths['star']/srr/'ReadsPerGene.out.tab', sep='\t', header=None, names=['Gene','Unstranded','Forward','Reverse'])
    g=g[~g['Gene'].str.startswith('N_',na=False)]; u=int(g['Unstranded'].sum()); f=int(g['Forward'].sum()); r=int(g['Reverse'].sum()); den=f+r
    ff=f/den if den else 0; rf=r/den if den else 0
    if rf>=0.80: call='REVERSE_STRANDED'; s=2
    elif ff>=0.80: call='FORWARD_STRANDED'; s=1
    elif max(ff,rf)<=0.65: call='UNSTRANDED'; s=0
    else: call='AMBIGUOUS_REVIEW'; s=-1
    strand_rows.append([srr,u,f,r,ff,rf,call,s])
mapq=pd.DataFrame(map_rows, columns=['SRR','Input_reads','Average_input_read_length','Uniquely_mapped_reads','Uniquely_mapped_pct','Average_mapped_length','Mismatch_rate_per_base_pct','Multi_loci_pct','Too_many_loci_pct','Unmapped_too_short_pct','Unmapped_other_pct']).merge(runmap[['GSM','SRR','Contrast_id','Role']],on='SRR',how='left').sort_values(['GSM','SRR'])
strand=pd.DataFrame(strand_rows, columns=['SRR','Gene_counts_unstranded','Gene_counts_forward','Gene_counts_reverse','Forward_fraction','Reverse_fraction','Strand_call','featureCounts_s']).merge(runmap[['GSM','SRR']],on='SRR',how='left').sort_values(['GSM','SRR'])
mapq.to_csv(os_paths['starqc']/f'{GSE}_STAR_mapping_QC_v1.tsv',sep='\t',index=False)
strand.to_csv(os_paths['starqc']/f'{GSE}_STAR_strandedness_v1.tsv',sep='\t',index=False)
if (strand['featureCounts_s']<0).any():
    bad = strand.loc[strand['featureCounts_s'] < 0].copy()
    bad.to_csv(os_paths['meta']/f'{GSE}_ambiguous_strandedness_affected_v1.tsv', sep='\t', index=False)
    processing_stop('REVIEW_REQUIRED_AMBIGUOUS_STRANDEDNESS', 23, [f'AFFECTED_SRR={len(bad)}'])
if strand.groupby('GSM')['featureCounts_s'].nunique().max()>1:
    discordant = strand.groupby('GSM')['featureCounts_s'].nunique()
    discordant = discordant[discordant > 1]
    pd.DataFrame({'GSM': discordant.index, 'Distinct_strand_calls': discordant.values}).to_csv(
        os_paths['meta']/f'{GSE}_within_GSM_strand_discordance_v1.tsv', sep='\t', index=False
    )
    processing_stop('REVIEW_REQUIRED_WITHIN_GSM_STRAND_DISCORDANCE', 24, [f'AFFECTED_GSM={len(discordant)}'])

# -------------------- featureCounts per SRR --------------------
fc_tables=[]; fc_qc=[]
input_read_map=mapq.set_index('SRR')['Input_reads'].to_dict()
for row in strand.itertuples(index=False):
    srr=row.SRR; sval=int(row.featureCounts_s); bam=os_paths['star']/srr/'Aligned.sortedByCoord.out.bam'; out=os_paths['fc']/f'{srr}_featureCounts_v1.txt'
    if not out.exists() or not Path(str(out)+'.summary').exists():
        cmd=[str(FC),'-T',str(args.fc_threads),'-p','--countReadPairs','-B','-C','-s',str(sval),'-t','exon','-g','gene_id','-a',str(GTF),'-o',str(out),str(bam)]
        with open(os_paths['log']/f'{srr}_featureCounts_v1.log','w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    t=pd.read_csv(out,sep='\t',comment='#'); countcol=[c for c in t.columns if c not in ['Geneid','Chr','Start','End','Strand','Length']][0]
    sub=t[['Geneid',countcol]].rename(columns={countcol:srr}); fc_tables.append(sub)
    assigned=int(sub[srr].sum()); summ=pd.read_csv(str(out)+'.summary',sep='\t',index_col=0); scol=summ.columns[0]
    fc_qc.append([srr,sval,assigned,int(input_read_map[srr]),assigned/int(input_read_map[srr]),int(summ[scol].sum()),int(summ.loc['Unassigned_MultiMapping',scol]),int(summ.loc['Unassigned_NoFeatures',scol]),int(summ.loc['Unassigned_Ambiguity',scol])])

srr_counts=fc_tables[0]
for t in fc_tables[1:]: srr_counts=srr_counts.merge(t,on='Geneid',how='outer',validate='one_to_one')
srr_counts=srr_counts.fillna(0)
for srr in expected_srr: srr_counts[srr]=srr_counts[srr].astype(np.int64)
assert srr_counts.shape[0]==78348
srr_counts=srr_counts[['Geneid']+expected_srr]
srr_counts.to_csv(os_paths['gsm']/f'{GSE}_SRR_counts_v1.tsv',sep='\t',index=False)
fcqc=pd.DataFrame(fc_qc,columns=['SRR','featureCounts_s','Assigned_fragments','STAR_input_read_pairs','Assigned_per_input_pair','FC_summary_total_alignments','FC_Unassigned_MultiMapping_alignments','FC_Unassigned_NoFeatures','FC_Unassigned_Ambiguity']).merge(runmap[['GSM','SRR','Contrast_id','Role']],on='SRR',how='left').sort_values(['GSM','SRR'])
fcqc.to_csv(os_paths['starqc']/f'{GSE}_featureCounts_QC_v1.tsv',sep='\t',index=False)

# -------------------- SRR -> GSM aggregation --------------------
gsm_counts=pd.DataFrame({'Geneid':srr_counts['Geneid']}); agg=[]; corr=[]; logc=np.log1p(srr_counts.set_index('Geneid'))
for gsm,z in runmap.groupby('GSM'):
    srrs=sorted(z['SRR'].tolist()); gsm_counts[gsm]=srr_counts[srrs].sum(axis=1).astype(np.int64)
    before=int(srr_counts[srrs].to_numpy().sum()); after=int(gsm_counts[gsm].sum()); role=z['Role'].unique()[0]; contrast=z['Contrast_id'].unique()[0]
    agg.append([gsm,role,contrast,len(srrs),','.join(srrs),before,after,after-before,'PASS' if after==before else 'FAIL'])
    if len(srrs)>1:
        c=logc[srrs].corr(); vals=[c.iloc[i,j] for i in range(len(srrs)) for j in range(i+1,len(srrs))]
        corr.append([gsm,len(srrs),min(vals),float(np.mean(vals)),max(vals)])
    else:
        corr.append([gsm,1,np.nan,np.nan,np.nan])
agg=pd.DataFrame(agg,columns=['GSM','Role','Contrast_id','SRR_count','SRRs','SRR_total_assigned','GSM_total_after_sum','Difference','Aggregation_QC']).sort_values('GSM')
if not (agg['Aggregation_QC']=='PASS').all():
    agg.to_csv(os_paths['starqc']/f'{GSE}_SRR_to_GSM_aggregation_QC_v1.tsv',sep='\t',index=False)
    processing_stop('SRR_TO_GSM_AGGREGATION_FAIL', 25, [f'FAILED_GSM={int((agg["Aggregation_QC"]!="PASS").sum())}'])
corr=pd.DataFrame(corr,columns=['GSM','SRR_count','Min_pairwise_log1p_Pearson','Mean_pairwise_log1p_Pearson','Max_pairwise_log1p_Pearson']).sort_values('GSM')
gsm_counts=gsm_counts[['Geneid']+expected_gsm]
gsm_counts.to_csv(FINAL,sep='\t',index=False)
agg.to_csv(os_paths['starqc']/f'{GSE}_SRR_to_GSM_aggregation_QC_v1.tsv',sep='\t',index=False)
corr.to_csv(os_paths['starqc']/f'{GSE}_within_GSM_SRR_concordance_v1.tsv',sep='\t',index=False)
runmap.sort_values(['Contrast_id','Role','GSM','SRR']).to_csv(os_paths['gsm']/f'{GSE}_sample_map_v1.tsv',sep='\t',index=False)

# Final checks and status
assert gsm_counts.shape[0]==78348 and len(gsm_counts.columns)-1==len(expected_gsm)
summary = [
    f'GSE={GSE}', f'GSM={len(expected_gsm)}', f'SRR={len(expected_srr)}', f'FASTQ={len(x)}',
    f'MIN_ACTUAL_READ_LENGTH={int(fq.Min_read_length.min())}', f'MAX_ACTUAL_READ_LENGTH={int(fq.Max_read_length.max())}',
    f'MIN_UNIQUE_MAPPING_PCT={mapq.Uniquely_mapped_pct.min():.4f}', f'MAX_UNIQUE_MAPPING_PCT={mapq.Uniquely_mapped_pct.max():.4f}',
    'STRAND_CALLS=' + ','.join(f'{k}:{v}' for k,v in strand.Strand_call.value_counts().to_dict().items()),
    f'GSM_AGGREGATION_PASS={int((agg.Aggregation_QC=="PASS").sum())}/{len(agg)}', 'STATUS=PASS'
]
DONE.write_text('\n'.join(summary)+'\n')
write_status(
    'PASS',
    [
        f'GSM={len(expected_gsm)}',
        f'SRR={len(expected_srr)}',
        f'FINAL_GSM_MATRIX={FINAL.relative_to(ROOT)}'
    ]
)

# checksums for final deliverables
files=[FINAL, os_paths['gsm']/f'{GSE}_sample_map_v1.tsv', os_paths['starqc']/f'{GSE}_STAR_mapping_QC_v1.tsv', os_paths['starqc']/f'{GSE}_STAR_strandedness_v1.tsv', os_paths['starqc']/f'{GSE}_featureCounts_QC_v1.tsv', os_paths['starqc']/f'{GSE}_SRR_to_GSM_aggregation_QC_v1.tsv', DONE]
with open(os_paths['meta']/f'{GSE}_final_outputs_SHA256_v1.txt','w') as o:
    for p in files:
        o.write(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(ROOT))+'\n')

if args.cleanup:
    deleted=[]
    for p in sorted(os_paths['raw'].glob('*.fastq.gz')):
        deleted.append([str(p.relative_to(ROOT)),p.stat().st_size,'FASTQ']); p.unlink()
    for p in sorted(os_paths['star'].glob('*/Aligned.sortedByCoord.out.bam')):
        deleted.append([str(p.relative_to(ROOT)),p.stat().st_size,'BAM']); p.unlink()
    pd.DataFrame(deleted,columns=['Deleted_path','Bytes','Type']).to_csv(os_paths['meta']/f'{GSE}_cleanup_manifest_v1.tsv',sep='\t',index=False)

print(f'===== {GSE} PRODUCTION PASS =====')
print('GSM matrix:', FINAL)
print('GSM:',len(expected_gsm),'SRR:',len(expected_srr))
print('Cleanup:', 'YES' if args.cleanup else 'NO')
