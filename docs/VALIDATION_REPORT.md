# Server-side validation evidence — 2026-10-13

The original project root was read-only and new outputs were generated in separate validation directories. The selected release scripts were invoked from the v0.2 code package; v0.3 keeps the same **89 source script bytes** and updates documentation only.

| Stage | Code/input validation | Result | Important limitation |
|---|---|---|---|
| 1 | 89 checksums, 66 Python, 5 shell, 18 R parse; R and Python packages | PASS (all syntax errors = 0) | STAR and featureCounts were not on PATH; external executables may require explicit paths |
| 2A | Stored scientific outputs vs expected manuscript numbers | 30 PASS; 0 failures | Stored-output audit, no model refit |
| 2B | R1C1 threshold, R1C5 heterogeneity and R1C6 confounding scripts | PASS | Reads frozen upstream model output; human cohort exclusion omitted from this specific run |
| 2C | 12 sampled mouse multilevel REML models (4 core, 4 significant noncore, 4 nonsignificant) | 156 compared numerical fields, 0 discrepancies | Subset only; no genome-wide FDR |
| 2D | `ENSMUSG00000129703.1` BFGS rescue refit | 13 numerical fields, 0 discrepancies | Subset only; no genome-wide FDR |
| 3A | Whole pooled Hallmark GSEA | 50 pathways; 25 FDR-significant; same significant pathway set | NES/p-value compared with numeric/stochastic tolerance |
| 3B | Whole duration and fat-adjusted duration Hallmark GSEA | 50 pathways in each; 29 and 25 significant; both lists matched | NES/p-value compared with numeric/stochastic tolerance |
| 3C | 5 sampled genes, duration/fat/control moderator, matched-null, joint models | 510 fields, 0 discrepancies | Subset only; no global BH-FDR |
| 3D | 12 sampled human genes, primary RE/FE and strict4 sensitivity | 456 fields, 0 discrepancies | Subset only; no genome-wide BH-FDR |
| 3E | Whole mouse–human ortholog joining, directional concordance and 8 logistic models | 10,658 orthologs, 5,484 same-direction, primary OR=1.074563390359; **426,487 compared non-key cells, 0 discrepancies** | Frozen upstream mouse/human meta results and ortholog map used; no upstream model refits in this stage |

Additional separate verification: R1C14 leave-one-human-cohort-out logistic FE analysis against the submitted workbook returned the expected 5 model rows. The 105-contrast apeglm sensitivity output was inspected (57/59 core genes retained the directional consistency cutoff), but this report does not claim an independent full apeglm production rerun.

### Explanation of the Stage 3E `numeric_parsed=0` log field

The comparison function checks `a == b` **before** trying `float(a)` and `float(b)`. If a numeric field has precisely the same textual representation in the stored and re-executed output, it is counted as a direct string equality, not as `numeric_parsed`. Equivalent missing-value markers also count as matching. Therefore the zero counter represents an unnecessary numeric-conversion fallback, **not** zero numeric fields included in comparison. All 426,487 non-key fields were visited, and zero discrepancies were reported.

### Boundaries of what was demonstrated

Re-running a *released analysis script* against frozen, locally available upstream datasets supports reproducibility of that analytical step. It is not a fresh end-to-end reproduction from GEO/SRA; original selection/screening adjudication and full-sequence pre-processing have not been independently re-run. Code and documentation should not imply otherwise. For public use, provide archive links or clear acquisition instructions for upstream frozen data, software versions, genome/annotation references, and licensed gene sets.
