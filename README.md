# HFD mouse liver RNA-seq meta-analysis — analysis code and reproducibility documentation

> This repository provides selected original analysis scripts and frozen metadata for a meta-analysis of 64 mouse liver GEO Series (105 HFD-versus-control contrasts) and three human liver cohorts. Selected model refits and complete downstream Hallmark/concordance analyses were verified against frozen inputs. **A fresh end-to-end run from raw FASTQ files has not been independently demonstrated.**

## Associated frozen data (Zenodo)

Frozen analysis inputs, derived outputs, sensitivity results, and reference annotations are deposited separately as four data archives in a Zenodo Dataset. **Reserved Dataset DOI:** `10.5281/zenodo.23355111`. After publication, the dataset can be accessed at `https://doi.org/10.5281/zenodo.23355111`.

- `HFD_01_analysis_inputs_v1.zip`: 15 frozen analytical input files, including full contrast-level meta-analysis effects.
- `HFD_02_analysis_results_v1.zip`: 33 analysis result files.
- `HFD_03_reviewer_sensitivity_v1.zip`: 17 sensitivity-analysis outputs.
- `HFD_04_reference_annotations_v1.zip`: frozen MSigDB mouse Hallmark gene sets and Ensembl release 116 ortholog mappings, with source/rights documentation.

The dataset also provides `ZENODO_ARCHIVE_SHA256.tsv` and `ZENODO_SOURCE_SHA256.tsv` for integrity checks. The code and data archives have complementary scopes; raw FASTQ files, full count matrices, and some intermediate files are not included. Raw reads remain available from the source GEO/SRA accessions.

## What is included

- `geo_search/`: historical GEO search terms and screening-support scripts. 
- `scripts/`: selected 107-candidate GEO sample–run mapping, acquisition, RNA-seq production, manual-metadata adjudication, mouse DESeq2/multilevel REML/core, moderator, Hallmark pathway, human meta-analysis and ortholog/concordance scripts. Original versioned filenames are preserved.
- `reviewer_sensitivity/`: **8** additional revision analyses, including threshold/LOO, heterogeneity, GSE-aggregate REML, moderator confounding, human cohort-exclusion FE logistic sensitivity, and apeglm shrinkage.
- `00_metadata/`: selected frozen accession, study-design, and curation metadata. Large frozen analytical inputs and result tables are provided separately through the Zenodo Dataset; complete raw/count and intermediate inputs are not bundled.
- `expected_results/`: previously verified small R1C14 cohort-exclusion output (for comparison).
- `docs/`: execution stages, dependencies, SHA-256 provenance, portability changes, and missing inputs.

**Figure-rendering scripts are intentionally excluded.** Scientific analysis and quantitative outputs are the scope; graphical reproduction of PDF figures is not claimed.

## Running from the correct project root

Create a project directory matching the source `00_metadata/`, `03_counts/`, `04_de/`, `05_meta/`, `06_pathway/`, and `07_human_validation/` folder structure. From that directory, set:

```bash
export HFD_PROJECT_ROOT="$(pwd)"
python3 tools/check_frozen_inputs.py --root "$HFD_PROJECT_ROOT"
```

Where STAR or featureCounts is used, supply absolute executable paths using `export HFD_STAR=/path/to/STAR` and `export HFD_FEATURECOUNTS=/path/to/featureCounts` (production script checks existence). Then run selected scripts in the scientific order shown in `docs/WORKFLOW.md`. This is **not** a turnkey `run-all` pipeline: raw sequencing reprocessing requires substantial storage and curated intermediate inputs. GEO search commands require NCBI EDirect; sequencing steps require STAR, featureCounts, and other installed software. See `docs/INPUT_REQUIREMENTS.md` and `docs/RELEASE_CHECKLIST.md`.

### Standalone verified reviewer sensitivity

With the exact submitted supplementary workbook (not bundled in this code repository):

```bash
python3 reviewer_sensitivity/verify_R1C14_human_cohort_LOO_sensitivity.py \
  --workbook /path/to/SupplementaryTables_v20260922.xlsx \
  --output outputs/R1C14_human_cohort_exclusion.csv
```

The primary comparison uses three-human-cohort **RE** estimates already present in the workbook. Each leave-one-human-cohort-out comparison recombines the **remaining two cohorts by inverse-variance fixed effects**, *not* two-cohort REML. This code only evaluates logistic directional concordance.

For the core59 `apeglm` sensitivity, run with the full-contrast argument `105`, not the script's two-contrast pilot default:

```bash
Rscript reviewer_sensitivity/run_R2C1_LFC_shrinkage_sensitivity_v1.R 105
```

## Supplementary Tables

`scripts/135_build_canonical_supplementary_tables_v1_1.py` generates historical expanded source Tables S1–S22. The submitted **S1–S17** workbook was subsequently **manually selected and assembled** from these source results. It is not claimed that a single script reproduces the final workbook layout exactly. See `docs/SUPPLEMENTARY_TABLES.md`.

## Caveats

- The frozen `HFD_primary_analysis_GSE_frozen_v1.txt` has **107** candidates, *not* the final 64 GSEs. The final 64 GSEs are defined by the frozen **105-contrast** audit/membership data and original study selection records.
- The code repository includes small metadata only. A separate Zenodo Dataset contains selected full gene-level analytical inputs/results and frozen MSigDB/Ensembl references. Raw FASTQ files, count matrices, full human source matrices and some intermediate files still require retrieval or regeneration.
- ROOT portability patches change path configuration only, not model formulas or eligibility logic. Selected code paths passed server-side re-execution checks, but portability and the entire raw-data workflow were not independently validated in a clean environment. `docs/PATH_PATCH_REPORT.tsv` records modified files; `docs/RELEASE_SHA256.tsv` records the curated source-code checksums.
- GEO discovery dates and candidate-set selection should not be inferred from rerunning live NCBI queries, whose results may change.
- The frozen reference datasets have independent provenance and reuse terms; see the source/rights note in the Zenodo reference archive. Do not assume all underlying third-party materials share the same license as this repository.

## Re-execution evidence and current status (2026-10-13)

- Server-side Stage 1: **89/89 code-file hashes matched**, Python 66/66, Bash 5/5, and R 18/18 syntax checks passed; principal model packages detected on the server.
- Stage 2A: **30/30 existing-output checks** agreed with report values. This is output verification, not model re-execution.
- Stage 2B: released R1C1 threshold, R1C5 heterogeneity, R1C6 metadata-confounding scripts successfully **re-executed** with expected numbers; R1C14 was not included in that run (it had been checked separately against the submitted workbook).
- Stages 2C–2D: **12 primary multilevel REML gene fits** and **one BFGS rescue** matched archived parameters (156 + 13 fields, 0 discrepancies). Genome-wide BH-FDR was deliberately excluded from this subset check.
- Stages 3A–3B: complete **pooled Hallmark (25 significant)**, **duration Hallmark (29)** and **fat-adjusted duration Hallmark (25)** results re-executed and significance lists agreed; numerical/stochastic tolerances were used for GSEA scores.
- Stage 3C: **five selected genes** refitted across duration/fat/control, corresponding matched-null and joint moderator models; **510 compared fields, 0 discrepancies**.
- Stage 3D: **12 human genes** refitted with the public code, including strict-control-4 sensitivity; **456 compared fields, 0 discrepancies**.
- Stage 3E: full downstream **10,658-ortholog concordance and eight logistic models** re-executed; 5,484 same-direction genes and primary OR 1.074563390359 reproduced. **426,487 non-key table cells** matched under the comparator (exact text match or equivalent missing representations, with numeric tolerance when required). Its `numeric_parsed=0` diagnostic means the direct-match shortcut resolved all comparable entries; it does **not** indicate that numbers were omitted from the comparison.

**These tests establish selected model and downstream analytical reproducibility against frozen inputs.** They do **not** establish the reproducibility of raw FASTQ acquisition, preprocessing, the full 105-contrast DESeq2 workflow, genome-wide REML refitting or independent reconstruction of manual study screening. The repository is an **analytical code release with explicitly limited reproducibility scope**, not a turnkey workflow or an independent FASTQ-to-manuscript reproduction. Some documents under `docs/` preserve historical pre-release audit snapshots; follow this README for current release scope.

See `docs/VALIDATION_REPORT.md` for stage-by-stage evidence. `docs/RELEASE_GAPS.md` and `docs/RELEASE_CHECKLIST.md` retain a historical audit of scope limitations and pre-release tasks.
