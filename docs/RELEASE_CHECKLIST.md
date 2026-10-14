# Publication gating checklist — updated 2026-10-13

## Confirmed so far

- [x] Curate 89 source files (73 `scripts/`, 7 GEO-search, 8 reviewer, 1 input tool) and preserve original filenames. No figure-rendering scripts.
- [x] Replace old absolute project root with `HFD_PROJECT_ROOT` on *release copies* of 67 source scripts; retain code checksums and patch report.
- [x] Pass checksums, Python/Bash/R syntax and basic server package checks (Stage 1).
- [x] Compare stored outputs to manuscript figures (30 cross-checks, Stage 2A).
- [x] Re-execute three reviewer-analysis scripts and match expected values (Stage 2B).
- [x] Refit selected mouse REML models, including one BFGS rescue (Stages 2C/2D).
- [x] Re-execute three Hallmark pathway analyses (Stages 3A/3B).
- [x] Refit selected moderator and human meta models (Stages 3C/3D).
- [x] Re-execute the entire downstream 10,658-ortholog/8-logistic-model analysis (Stage 3E).
- [x] Verify R1C14 FE human-cohort-exclusion logistic script against saved output, independent of Stage 2B.
- [x] Explicitly document manually assembled submitted supplementary workbook.

## Must finish before public release

- [ ] Check every candidate script's input dependencies and actual downstream output link; document final execution commands, not just stage ordering.
- [ ] Freeze actual R/Python/system dependency versions and applicable reference data.
- [ ] Decide how the needed 106,366,306-byte meta-effect input and other frozen datasets will be accessed (data repository / rebuild instructions).
- [ ] Provide a complete final 64-GSE / 105-contrast selection and moderator design trail, with 671→107→64 provenance limits clearly stated.
- [ ] Document human phenotype/group assignments and cohort-level source matrix retrieval.
- [ ] Crosswalk the original generated S1–S22 tables to final S1–S17 by content (no scripted final layout required).
- [ ] Check code and redistributable metadata for credentials, privacy, personal paths and third-party licensing restrictions.
- [ ] Add a suitable license, citation metadata and dated release notes.
- [ ] Publish GitHub release and Zenodo DOI, then insert actual identifiers in manuscript/Data Availability/reviewer responses.

## Scope limitation

No independently re-executed full 105-contrast DESeq2 processing from raw FASTQ, nor full genome-wide REML fit of all genes, has been demonstrated. Stage 3E whole-output reproduction depended on frozen upstream results. Do not claim an end-to-end fresh reproduction.
