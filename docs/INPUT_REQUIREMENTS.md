# Required data and frozen metadata (selected frozen files now bundled)

| Input | Example project-relative source | Provenance/action |
|---|---|---|
| Historical screening-candidate GSE list | `00_metadata/frozen_search/HFD_primary_analysis_GSE_frozen_v1.txt` | **107 candidates**, not the final 64. The frozen 105-contrast analysis input identifies the final analysis universe. |
| FASTQ sample manifest | `00_metadata/technical_harmonization/frozen_v1/HFD_primary_ENA_R1R2_download_manifest_frozen_v1.tsv` | 951 SRA runs/837 GSM sample mapping; review redistributability |
| Final production QC | `00_metadata/processing_production/final_audit_v1/HFD_primary64_final_resolved_status_v1.tsv` | Required by `70` |
| Final moderator curation | `00_metadata/contrasts/moderators/manual_adjudication/HFD_primary105_moderator_manual_adjudication_batch11_v1.tsv` | Frozen human-curated metadata; preserve original evidence |
| Contrasts and memberships | `00_metadata/analysis_inputs/frozen_v1/*_v1_1.tsv` | 105 contrast memberships / annotations |
| Count matrix, read QC | `03_counts/`, `04_de/` and raw GEO/SRA FASTQ records | Rebuild from public sequencing data or provide accession-level acquisition instructions |
| Primary effect data | `00_metadata/analysis_inputs/meta_analysis_v1/HFD_primary105_meta_effects_primaryQC_v1.tsv.gz` | Gene-contrast source for primary model/reviewer analyses; large, not bundled |
| Primary gene meta output | `05_meta/gene_level_v1_2/HFD_gene_multilevel_meta_primary_v1_2.tsv.gz` | Source for thresholds, heterogeneity and sensitivity |
| Mouse Hallmark frozen set | `00_metadata/analysis_inputs/pathway_analysis_v1_2/MSigDB_mouse_Hallmark_frozen_v1_2.tsv.gz` | MSigDB licensing/redistribution must be checked; prefer versioned acquisition instructions |
| Human cohorts | `GSE126848`, `GSE130970`, `GSE135251`, per-sample group mapping and count matrices | Public GEO acquisition, cohorts processed with edgeR/limma |
| Orthology | Ensembl Compara release 116 high-confidence one-to-one mapping | Freeze source/version and mapping access instructions |
| Original supplementary source data | `09_supplementary/tables_v1_1/` | Earlier S1–S22; ensure submitted S1–S17 content crosswalk |

**Important:** `HFD_primary105_meta_effects_primaryQC_v1.tsv.gz` is about 106 MB by listed file size and should not be committed to a regular repository without checking size limits / a suitable data archive. The upstream FASTQ and count matrices should not be blindly committed either.

The bundled `00_metadata/` files provide a limited frozen snapshot and do not contain full raw FASTQ, final 64-study count matrices, human expression/count matrices, full model results, or external licensed gene sets. The source manuscript supplementary workbook is likewise not bundled.
