# Software, versions and reference sources

The release candidate does not contain a lockfile or container image. The following **versions/references are indicated by original scripts/manuscript**; confirm precise versions and commands on the actual analysis server and freeze a `sessionInfo()` / environment manifest before publication.

| Component | Source indication | Status |
|---|---|---|
| Mouse reference genome | GRCm39; GENCODE M39 GTF | Author's reference FASTA/GTF and STAR index not bundled |
| STAR | Original production executable path labelled 2.7.10b | Set `HFD_STAR=/absolute/path/to/STAR`; verify actual binary version |
| Subread featureCounts | Original production executable path labelled subread 2.0.6 | Set `HFD_FEATURECOUNTS=/absolute/path/to/featureCounts`; verify binary version |
| R | Manuscript provenance includes R 4.3.1 | Verify server `R --version` and `sessionInfo()` |
| DESeq2 | Version was recorded dynamically in original analysis output | Confirm installed version from production logs |
| metafor / data.table | Version recorded dynamically by scripts | Confirm installed versions from production logs |
| fgsea | Manuscript methods describe 1.28.0 | Confirm from package output; mouse MSigDB 2025.1.Mm Hallmark set must be obtained/licensed separately |
| edgeR / limma | Manuscript methods describe 4.0.16 / 3.58.1 | Confirm from cohort DE log |
| Orthology | Ensembl Compara release 116, high-confidence 1:1 | Freeze mapping and archive source provenance |
| Other R/Python libraries | Refer to `requireNamespace(...)` and imports in each selected script | Freeze a complete pip/conda and R package manifest |

## Pre-publication server commands

```bash
Rscript -e 'writeLines(capture.output(sessionInfo()), "R_sessionInfo.txt")'
STAR --version       # or "$HFD_STAR" --version
featureCounts -v     # or "$HFD_FEATURECOUNTS" -v
python3 --version
python3 -m pip freeze > python_requirements_frozen.txt
```

The `STAR` and `featureCounts` commands in this table are illustrative: use the actual installation paths and flags supported by installed binaries.
