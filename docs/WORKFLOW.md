# Analytical workflow — ordered by scientific dependency

**Not a turnkey executable shell script.** These are original selected filenames, not a validated run on portable paths. Within an ordered stage some scripts have additional data dependencies or require manual evidence adjudication.

| Stage | Relevant scripts / artifacts | Main analytical output |
|---|---|---|
| GEO discovery | `geo_search/00_discover_geo_hfd_liver.sh`, `geo_search/01_discover_geo_hfd_liver_strict_v2.sh`, `geo_search/02_rank_geo_candidates_v2.py`; retained original strict v1 and audit-pack scripts | Broad and strict query lists / triage only; curate frozen cohort separately |
| Acquisition and processing | `01_build_primary_sample_run_manifest.py` (107 candidates), `02_reconcile`, `03`, `04`, `06`, `07`, `20`, `29`, `30_*_v1_3_2.py`, `31_*_v1_4.sh`, `32`; resolve `66`–`69` QC exceptions | 64-GSE per-sample read counts / final QC status |
| Metadata adjudication | `40_*_v1_1_1.py`, `41`–`65` in substantive batch order; frozen batch11 human-reviewed annotations | Final 105-contrast moderator master |
| Mouse analysis inputs | `70_*_v1_1.py` → `71`, `72` | Frozen sample-to-contrast metadata and modelability |
| Mouse DESeq2 | `73_run_HFD_primary105_contrast_DESeq2_v1.R` | 105 contrast-level gene log2FC/SE |
| Multilevel meta | `74_build_HFD_primary105_meta_input_v1.py` → `75_run_HFD_gene_multilevel_meta_v1_1.R` → `76` audit → `77_retry_and_finalize_HFD_gene_multilevel_meta_v1_2.R` | Pooled gene-level REML effects, τ²/I², PI |
| Reproducible core | `78_finalize_HFD_stable_core_LOO_v1.R` → `79` annotation → `80` audit | 59-gene study-defined core |
| Mouse pooled pathway | `81_*_v1_2.R` → `82_*_v1_2.R` | Hallmark GSEA on complete gene ranking |
| Moderators and duration | `83_*_v1_2.R` → `84_*_v1.R` → `85`, `86`, `87` | Duration, fat, control-diet, integrated pathway estimates |
| Human cohorts | `88`, `89`, `90_*_v1_2.py` → `91_*_v1_2.R` → `92_*_v1.R` | Human cohort differential expression, 3-cohort REML/FE |
| Ortholog and concordance | `93_*_v1_6.R` → `94_*_v1_1.R` → `95_*_v1_1.R` → `96` | One-to-one mouse–human concordance models/pathways |
| Audits and provenance | `97`, `98`; `134` | Manuscript and supplementary-source consistency |
| Supplementary source tables | `135_build_canonical_supplementary_tables_v1_1.py` | Expanded analysis tables S1–S22 (NOT the assembled submitted S1–S17 workbook) |
| Reviewer results | `reviewer_sensitivity/IJMS_R1C1_*`, `IJMS_R1C5_*`, `IJMS_R1C6_*`, `run_R2C1_LFC_shrinkage_sensitivity_v1.R 105`, `verify_R1C14_human_cohort_LOO_sensitivity.py --workbook ...` | Sensitivity summaries used in reviewer responses |

The production wrapper explicitly invokes `30_process_one_GSE_production_v1_3_2.py`. Multilevel `77` refits are a continuation of `75` rather than an independent substitute.
